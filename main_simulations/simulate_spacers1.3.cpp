/* 21/04/25
Simple model for the dynamics of spacers in CRISPR arrays. One spacer is lost and another spacer is gained at each time step. The spacer that is lost is chosen with probability proportional to exp(-incidence), where the incidence (W) of the target depends on the time since the spacer acquisition. We asume that the incidence decays exponentially with time as

incidence: W(t) = log(1/a) * b^{-g} * b^t; with 0 < b < 1; 0 < a < 1; and g > 0. (It's ln the log)

This can also be written as W(t) = W0 * exp{-B * t}; with B = -log(b) and W0 = -log(a) * exp{B*g}. (W0 is the incidence at time 0.)

Assuming a Poisson process for spacer-target matching, the probability that a spacer matches its target per unit of time is Pmatch = 1 - exp(-W). The probability that there is no match (therefore the spacer is dispensable) is Pdisp = 1 - Pmatch = exp(-W). We use Pdisp as a weight to choose which spacer to remove.

Note that -log(b) determines the rate at which the incidence decays. The values of b are bounded between 0 and 1. Epidemic viruses will have small "b" (large B) and persistent viruses will have "b" close to 1 (B close to zero).

The prefactor b^{-g} generates a tradeoff between the initial incidence and the long-term incidence. The interpretation is that endemic (~persistent) viruses (b->1) produce sustained infections with low intensity, whereas epidemic (~lytic) viruses (b->0) lead to acute infections, with high but short-lived viral titers. Parameter "g" takes values from 0 to infinity (although it probably does not make sense to try values beyond 2 or 3). The greater the value of "g", the larger the difference between the initial incidence of epidemic and persistent viruses. If g = 0, epidemic and persistent viruses have the same initial incidences. In terms of the model, parameter "g" determines the time at which the incidence and removal probiabilities of all spacers coincide, regardless of "b". At t < g, selection favors spacers against epidemic viruses; at t > g, selection favors spacers against endemic viruses.

Parameter "a" takes values between 0 and 1 and determines the typical magnitude of viral incidence, with low "a" implying high incidences and high "a" implying low incidences. By tuning incidence values, "a" affects the probability that a target is matched by a cognate spacer. That becomes especially clear in the case of pure endemic targets (beta -> 1), whose probability of being matched is simply 1-a. Note that the removal probability for spacers against endemic viruses (before normalization) is equal to "a", whereas the long-term removal probability for spacers against epidemic viruses (before normalization) is 1. Therefore, in practice, parameter "a" determines the relative rates of removal of spacers against endemic and epidemic viruses in the long term. The lower the value of "a", the more favored the spacers against endemic viruses in the long term.

In this version of the model, parameters "a" and "g" are externally provided as an input. Parameter "b" is randomly chosen for each new spacer from a mixture probability distribution, such that with probability Q the paremeter "b" is set to 0.9999 (persistent virus), and with probability 1-Q, b is chosen from a uniform distribution (mixture of epidemic viruses with different retention levels). Thus, p(b) = 1-Q for b <= 0.999; p(0.999) = Q. The relative frequency of endemic targets (Q) is provided as external input.

The output of the simulation is Pmatch, the probability that the target of a spacer is present in the environment as a function of the spacer position in the array. The average values of "b", spacer age, and target incidence as a function of the position are also provided. Note that in the code we use the following nomenclature for these parameters:
	a -> alpha
	b -> beta
	g -> gamma
	Q -> pEndemic

For extreme parameter combinations (e.g. g=60 and Nsp=40), efficacy values could become very large, causing numerical overflow. 
To prevent this, in this script efficacy values were stored in log scale, and the fitness function is computed as  

𝑓𝑖𝑡𝑛𝑒𝑠𝑠=1−𝑒^−𝑒^log(𝑊)

only when efficacy is below a threshold (log(W) < 40). 
Above that, fitness is set to 1 directly, since the exponential of very large values is functionally indistinguishable from 1 in this context. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


#define MAXGEN 5000	// Number of generations to simulate
#define NUMSPC 40			// Number of spacers in the array
#define NUMREPS 5000		// Number of replicates

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;   // Required for function 'select()' which is used to obtain the quartiles.


//long seed = -1395170709;
long seed = - (long)time(NULL);	// Random seed for random number generation
long *idum = &seed;
double ran1(long *idum);

// Main functions
int SelectSpacer(double spEff[],int Nsp);						// Select a spacer with prob. proportional to its inefficacy
double CalcEfficacySingle(double logalpha,double beta,double gamma,double t);	// Decay in efficacy as a function of age
void UpdateArray(double spBeta[],int spAge[],int id,int Nsp,double pEndemic);	// Update array by replacing selected spacer
void UpdateEfficacy(double spBeta[],int spAge[],double spEff[],double logalpha,double gamma,int Nsp);	// Update the efficacies
void QuartilesByColumn(double Qs[][4],double data[][NUMSPC],int nrow,int ncol);				// Calculates quartiles (Q1, Q2, Q3) and mean
double ran1(long *idum);
double select(unsigned long k, unsigned long n, double arr[]);		// Used to calculate the quartiles
void init_array(double arr[],double val,int nrow);
void CalculateAveragesByColumn(double endemicAvgBySpacer[NUMSPC], double epidemicAvgBySpacer[NUMSPC],double data[][NUMSPC],int nrow,int ncol);	; //Calculate beta average for endemic and epidemic in a position of the spacer

int main (int argc,char *argv[])
{
  if (argc<4)
    {
    printf("Error: provide parameters: alpha, gamma, pEndemic\n");
    return -1;
    }


   // Define variables
   double AllBeta[NUMREPS][NUMSPC] = {0.};	// Betas for the spacers in all replicates
   int AllAge[NUMREPS][NUMSPC] = {0}; 	// Age of spacers in all replicates
   double AllEff[NUMREPS][NUMSPC] = {0.};	// Efficacy of spacers in all replicates
   double AllFitness[NUMREPS][NUMSPC] = {0.};	// Fitness of spacers in all replicates (1 - exp(efficacy))
   double QuartilesBeta[NUMSPC][4] = {0.}; 	// Quartiles (Q1, Q2, Q3) and mean beta for each position
   double QuartilesEff[NUMSPC][4] = {0.}; 	// Quartiles (Q1, Q2, Q3) and mean efficacy for each position
   double QuartilesFitness[NUMSPC][4] = {0.}; // Quartiles (Q1, Q2, Q3) and mean fitness for each position
   double QuartilesAge[NUMSPC][4] = {0.}; 	// Quartiles (Q1, Q2, Q3) and mean age for each position
   double TmpAge[NUMREPS][NUMSPC] = {0.};
   double endemicAvgBySpacer[NUMSPC] = {0.0};
   double epidemicAvgBySpacer[NUMSPC] = {0.0};
   double alpha, logalpha;		  // Probability of encounter with persistent virus
   double gamma = 5.;		  // Critical time at which all spacers are equally efficient regardless of beta
   double pEndemic = 0.2;	  // Relative frequency of endemic parasites
   int Nsp = NUMSPC;		  // Number of spacers in the array
   int Tmax = MAXGEN;		  // Number of time steps
   int Nrep = NUMREPS;		  // Number of replicates
   int ii,ii2;
   int validSimulation = 0;      // Flag to determine if the simulation is valid


   // Define functions
   void SingleSimulation(double spBetaOut[],int spAgeOut[],double logalpha,double gamma,double pEndemic,int Nsp,int tmax);

   // Parameter values: read from input
   alpha = atof(argv[1]);
   logalpha = -log(alpha);
   gamma = atof(argv[2]);
   pEndemic = atof(argv[3]);
    while (!validSimulation)
    {
        // Run simulations
        for (ii = 0; ii < Nrep; ii++)
            SingleSimulation(AllBeta[ii], AllAge[ii], logalpha, gamma, pEndemic, Nsp, Tmax);

        // Calculate final efficacies
        for (ii = 0; ii < Nrep; ii++)
            UpdateEfficacy(AllBeta[ii], AllAge[ii], AllEff[ii], logalpha, gamma, Nsp);

        // Calculate short-term final fitness
        for (ii = 0; ii < Nrep; ii++)
        {
            for (ii2 = 0; ii2 < Nsp; ii2++)
               if (AllEff[ii][ii2] < 40) // Avoid computing exp(-exp(AllEff)) for large values
                    AllFitness[ii][ii2] = 1. - exp(-exp(AllEff[ii][ii2])); // Reverse the logarithm when possible
               else
                    AllFitness[ii][ii2] = 1; // If spEff is too large, thhe fitness value is 1
        }

        // Calculate statistics
        QuartilesByColumn(QuartilesEff, AllEff, Nrep, Nsp);      // Quartiles and mean for efficacy
        QuartilesByColumn(QuartilesBeta, AllBeta, Nrep, Nsp);    // Quartiles and mean for beta
        QuartilesByColumn(QuartilesFitness, AllFitness, Nrep, Nsp); // Quartiles and mean for fitness
        for (ii = 0; ii < Nrep; ii++)
        {
            for (ii2 = 0; ii2 < Nsp; ii2++)
                TmpAge[ii][ii2] = AllAge[ii][ii2];
        }
        QuartilesByColumn(QuartilesAge, TmpAge, Nrep, Nsp);      // Quartiles and mean for age
        // Calculate beta average for endemic vs epidemic
        CalculateAveragesByColumn(endemicAvgBySpacer, epidemicAvgBySpacer, AllBeta, Nrep, Nsp);

        // Verify condition: last spacer's age < MAXGEN / 5 (in each replicate)
        int maxAge = 0;
        for (ii = 0; ii < Nrep; ii++)
        {
            if (AllAge[ii][Nsp - 1] > maxAge)
                maxAge = AllAge[ii][Nsp - 1];
        }
        if (maxAge < Tmax / 5)
        {
            validSimulation = 1; // Simulation is valid
        }
        else
        {
            Tmax *= 2; // Increase Tmax and repeat
        }
    }

    // Print results
    printf("alpha = %6.4f\n", alpha);
    printf("gamma = %6.4f\n", gamma);
    printf("prob.endemic = %6.4f\n", pEndemic);
    printf("Num. spacers = %i, Tmax = %i, Num. replicates = %i\n", Nsp, Tmax, Nrep);
    printf("seed = %li\n", seed);
	printf("beta\tage\tefficacy\tfitness\tendemicBeta\tepidemicBeta (means)\n");
	for (ii = 0; ii < Nsp; ii++)
    	printf("%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n", 
           QuartilesBeta[ii][3], 
           QuartilesAge[ii][3], 
           QuartilesEff[ii][1], 
           QuartilesFitness[ii][3], 
           endemicAvgBySpacer[ii], 
           epidemicAvgBySpacer[ii]);
    return 0;
}

// Auxiliary functions

void SingleSimulation(double spBetaOut[],int spAgeOut[],double logalpha,double gamma,double pEndemic,int Nsp,int tmax)
{
// Define variables
double spBeta[Nsp];	 // Betas for the spacers (perhaps initialize randomly)
int spAge[Nsp]; 	 	 // Age of spacers (perhaps initialize randomly)
double spEff[Nsp];	 // Efficacy of spacers
int idtoreplace = 0;	  	 // Index of the spacer that will be replaced
int tt = 0;		  	 // Time step

// Initialize arrays
init_array(spBeta,0.5,Nsp);
memset( spAge, 0, sizeof(spAge) );
// Initialize spAge with zeros that way, because I get an error in my compiler if I do int spAge[numSp]={0};

// Spacer dynamics
for (tt=0;tt<tmax;tt++)
   {
   UpdateEfficacy(spBeta,spAge,spEff,logalpha,gamma,Nsp);	// Update efficacies
   idtoreplace = SelectSpacer(spEff,Nsp);					// Select spacer to replace
   UpdateArray(spBeta,spAge,idtoreplace,Nsp,pEndemic);		// Update array
   }

// Assign to output vector
for (tt=0;tt<Nsp;tt++)
   {
   spBetaOut[tt] = spBeta[tt];
   spAgeOut[tt] = spAge[tt];
   }

return;
}

void init_array(double arr[],double val,int nrow)
{
int ii;
for (ii=0;ii<nrow;ii++)
   arr[ii] = val;
return;
}

void UpdateEfficacy(double spBeta[],int spAge[],double spEff[],double logalpha,double gamma,int Nsp)
{
int ii = 0;
for (ii=0;ii<Nsp;ii++)
   spEff[ii] = CalcEfficacySingle(logalpha,spBeta[ii],gamma,spAge[ii]);
return;
}

double CalcEfficacySingle(double logalpha,double beta,double gamma,double t)
{
double eff = log10(logalpha) + (-gamma * log10(beta)) + (t * log10(beta)); // it really is log(eff), page 122 of the notebook
return eff;
}

int SelectSpacer(double spEff[], int Nsp) {
    int ii = 0;
    double sumEff = 0.;
    double rnum;

    // Calculate total replacement susceptibility
    for (ii = 0; ii < Nsp; ii++) {
        if (spEff[ii] < 40) {  // Avoid computing exp(-exp(spEff)) for large values
            sumEff += exp(-exp(spEff[ii]));
        } else {
            sumEff += 0;  // If spEff is too large, its contribution is 0
        }
    }

    // Generate a random number in the range [0, sumEff)
    do {
        rnum = ran1(idum);
    } while (rnum >= 1);
    rnum *= sumEff;

    // Select spacer based on exp(-exp(spEff[ii]))
    sumEff = 0;
    ii = 0;
    do {
        if (spEff[ii] < 40) {  // Avoid unnecessary calculations
            sumEff += exp(-exp(spEff[ii]));
        }
        ii++;
    } while (sumEff < rnum && ii < Nsp);
    ii--;

    return ii;
}

void UpdateArray(double spBeta[],int spAge[],int id,int Nsp,double pEndemic)
{
int ii;
double getNewBeta(double pEndemic);

// Increase age of spacers after the removed one
for (ii=id+1;ii<Nsp;ii++)
   spAge[ii] = spAge[ii]+1;

// Replace spacers from the removed one to the first one by the previous spacer and increase age
for (ii=id-1;ii>=0;ii--)
   {
   spBeta[ii+1] = spBeta[ii];
   spAge[ii+1] = spAge[ii]+1;
   }

// Add new spacer
spAge[0] = 0;
spBeta[0] = getNewBeta(pEndemic);
return;
}

double getNewBeta(double pEndemic)
{
double rnum;

do rnum  = ran1(idum);
while (rnum >= 1);

// Biased bimodal distribution: epidemic vs endemic viruses
if (rnum >= pEndemic)
   return (rnum-pEndemic)/(1.-pEndemic);
else
   return 0.9999;

}

void QuartilesByColumn(double Qs[][4],double data[][NUMSPC],int nrow,int ncol)
{
double datacol[nrow];
double rowsum = 0;
int irow = 0;
int icol = 0;

for (icol=0;icol<ncol;icol++)
   {
   rowsum = 0;
   for (irow=0;irow<nrow;irow++)
      {
      rowsum += data[irow][icol];
      datacol[irow] = data[irow][icol];		// Populate datacol
      }
   Qs[icol][0] = select(nrow/4 - 1,nrow,datacol);	// Calculate quartiles
   Qs[icol][1] = select(nrow/2 - 1,nrow,datacol);
   Qs[icol][2] = select(nrow - nrow/4 - 1,nrow,datacol);
   Qs[icol][3] = rowsum / nrow;			// Calculate mean
   }

return;
}

void CalculateAveragesByColumn(double endemicAvgBySpacer[], double epidemicAvgBySpacer[], double data[][NUMSPC], int nrow, int ncol) {
double datacol[nrow];
double rowsum = 0;
int irow = 0;
int icol = 0;

for (icol=0;icol<ncol;icol++) {
        double endemicSum = 0.0, epidemicSum = 0.0;
        int endemicCount = 0, epidemicCount = 0;

        for (irow=0;irow<nrow;irow++) {
            if (fabs(data[irow][icol] - 0.9999) < 1e-4) {  // Endemic, take into account that C++ keep values like 0.9999 in  aproximate values
                endemicSum += data[irow][icol];
                endemicCount++;
            } else {  // Epidemic
                epidemicSum += data[irow][icol];
                epidemicCount++;
            }
        }
        // Beta epidemic and beta endemic for each position
        endemicAvgBySpacer[icol] = (endemicCount > 0) ? (endemicSum / endemicCount) : 0.0; //conditional to avoid dividing by zero
        epidemicAvgBySpacer[icol] = (epidemicCount > 0) ? (epidemicSum / epidemicCount) : 0.0;
    }
}

double select(unsigned long k, unsigned long n, double arr[])
/* Returns the kth smallest value in the array arr[1..n]. The input array will be rearranged
to have this value in location arr[k], with all smaller elements moved to arr[1..k-1] (in
arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).
*/
{
   unsigned long i,ir,j,l,mid;
   double a,temp;
   l=1;
   ir=n;
   for (;;) {
      if (ir <= l+1) { 			// Active partition contains 1 or 2 elements.
         if (ir == l+1 && arr[ir] < arr[l]) { // Case of 2 elements.
            SWAP(arr[l],arr[ir])
         }
         return arr[k];
      } else {
         mid=(l+ir) >> 1;			// Choose median of left, center, and right elements as partitioning element a. Also rearrange so that arr[l] ≤ arr[l+1], arr[ir] ≥ arr[l+1].
         SWAP(arr[mid],arr[l+1])
         if (arr[l] > arr[ir]) {
            SWAP(arr[l],arr[ir])
         }
         if (arr[l+1] > arr[ir]) {
            SWAP(arr[l+1],arr[ir])
         }
         if (arr[l] > arr[l+1]) {
            SWAP(arr[l],arr[l+1])
         }
         i=l+1;				// Initialize pointers for partitioning.
         j=ir;
         a=arr[l+1]; 				// Partitioning element.
         for (;;) {				// Beginning of innermost loop.
            do i++; while (arr[i] < a); 	// Scan up to find element > a.
            do j--; while (arr[j] > a);	// Scan down to find element < a.
            if (j < i) break;			// Pointers crossed. Partitioning complete.
            SWAP(arr[i],arr[j])
         }					// End of innermost loop.
         arr[l+1]=arr[j];			// Insert partitioning element.
         arr[j]=a;
         if (j >= k) ir=j-1;			// Keep active the partition that contains the kth element.
         if (j <= k) l=i;
      }
   }
}

double ran1(long *idum)
{
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   double temp;

   if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7; j>=0; j--) {
         k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR+k;
         if (*idum < 0) *idum += IM;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j] = *idum;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
