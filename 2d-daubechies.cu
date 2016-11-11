#include<stdio.h>
#include<math.h>
#include<sys/time.h>
#include<sys/time.h>
#include<stdlib.h>
#include"ext.h"
#define a 16
#define b 16

struct timeval startwtime, endwtime;
struct timeval startwtime1, endwtime1;
double seq_time;
int N,m,n;
void decompose(float **target, float *filter, float *quadfilter, float **input, int level);

__global__ void kernel1(float *target, float *input, float *filter, float *quadfilter, int m, int n, int megethos)  {
	
	__shared__ float imagesh[532];
	__shared__ float filt[20];
	int i;
	int starttemp=blockIdx.x*512;  //8
	int start=(starttemp/n)*(n+2*(megethos/2-1))+starttemp%n;
	int tid=threadIdx.x*blockDim.y+threadIdx.y;
	
	//metafora twn kymatidiwn stin shared memory
	if (tid<megethos) {
		if (blockIdx.y==1) {
			filt[tid]=filter[tid];
		}
		else {
			filt[tid]=quadfilter[tid];
		}
	}

	//metafora toy kommatioy toy pinaka pou tha xrisimopoiithei ap to sygkekrimeno block sthn shared memory
	imagesh[tid]=input[start+tid];
	imagesh[tid+256]=input[start+tid+256];
	if (tid+3<=megethos) {
		imagesh[tid+512]=input[start+tid+512];
	}
	

	__syncthreads();

	//ypologismoi
	float sum=0.0;
	for (i=0;i<megethos;i++) {
		sum+=filt[megethos-i-1]*imagesh[2*tid+i];
	}
	
	//eggrafi tou apotelesmatos
	int or1=(starttemp/2+tid)/(n/2);
	int or2=(starttemp/2+tid)%(n/2);
	if (blockIdx.y==1) {
		target[(or1+megethos/2-1)*n+or2]=sum;
	}else {
		target[(or1+megethos/2-1)*n+or2+n/2]=sum;
		
	}
}

__global__ void kernel2(float *target, float *target2, float *input, float *filter, float *quadfilter, int m, int n, int total2, int megethos)  {
	
	__shared__ float imagesh[532];
	__shared__ float filt[20];
	int or1, or2, i;
	int starttemp=blockIdx.x*512;
	int start=(starttemp/m)*n+starttemp%m;
	int tid=threadIdx.x*blockDim.y+threadIdx.y;
	
	//metafora twn kymatidiwn stin shared memory
	if (tid<megethos) {
		if (blockIdx.y==1 || blockIdx.y==3) {
			filt[tid]=filter[tid];
		}
		else {
			filt[tid]=quadfilter[tid];
		}
	}
	
	//metafora toy kommatioy toy pinaka pou tha xrisimopoiithei ap to sygkekrimeno block sthn shared memory
	if (blockIdx.y==1 || blockIdx.y==2) {
		or1=(start+tid)/n;
		or2=(start+tid)%n;
		imagesh[tid]=input[or2*n+or1];
		
		or2+=256;
		imagesh[tid+256]=input[or2*n+or1];
		if (tid+3<=megethos) {
			
			or2+=256;

			imagesh[tid+512]=input[or2*n+or1];
		}
	}else {

		or1=(start+tid)/n;
		or2=(start+tid)%n;
		imagesh[tid]=input[or2*n+or1+n/2];
		
		or2+=256;
		imagesh[tid+256]=input[or2*n+or1+n/2];
		if (tid+3<=megethos) {
			
			or2+=256;

			imagesh[tid+512]=input[or2*n+or1+n/2];
		}
	}

	__syncthreads();

	//ypologismoi
	float sum=0.0;
	for (i=0;i<megethos;i++) {
		sum+=filt[megethos-i-1]*imagesh[2*tid+i];
	}
	
	//eggrafi tou apotelesmatos
	or1=(starttemp/2+tid)/(m/2);
	or2=(starttemp/2+tid)%(m/2);
	if (blockIdx.y==1) {
		target[or2*total2+or1]=sum;
		target2[or2*(n/2+2*(megethos/2-1))+or1+megethos/2-1]=sum;
	}else if (blockIdx.y==2)
		target[or2*total2+or1+n/2]=sum;
	else if (blockIdx.y==3)
		target[(or2+m/2)*total2+or1]=sum;
	else
		target[(or2+m/2)*total2+or1+n/2]=sum;
	
}

		


int main(int argc,char **argv) {
	
	//arxikopoiisi kai memory allocation twn metavlhtwn pou tha xrisimopoiithoun
	float **imageserial, **imageoutserial, **usingimageserial;
	float **err;
	
	int k, sum1, sum2, sum3, sum4;
	float *image, *imageout, *usingimage, *devtempimage, *devimageinterout, *imageinterout, *devfilter, *devquadfilter, *devimageout;
	float *image1, *image2, *image3;
	int i,level, j, griddim1;
	float *filter, *quadfilter;
	level=atoi(argv[3]);
	N=atoi(argv[1]);
	m=atoi(argv[2]);
	n=m;
	filter=(float *)malloc(2*N*sizeof(float));
	quadfilter=(float *)malloc(2*N*sizeof(float));
	for (i=0;i<2*N;i++) {
		
		filter[i]=daubechies[N-1][i];
		quadfilter[i]=pow(-1.0,i)*daubechies[N-1][2*N-1-i];

	}
	imageserial=(float **)malloc(m*sizeof(float *));
	imageoutserial=(float **)malloc(m*sizeof(float *));
	for (i=0;i<m;i++)  {
		imageserial[i]=(float *)malloc(n*sizeof(float));
		imageoutserial[i]=(float *)malloc(n*sizeof(float));
	}
	image=(float *)malloc(m*n*sizeof(float));
	imageout=(float *)malloc(m*n*sizeof(float));
	
	for (i=0;i<m;i++)  {
		for (j=0;j<n;j++)  {
			image[i*n+j]=(float) rand()/RAND_MAX;
			imageserial[i][j]=image[i*n+j];
		}
	}

	usingimageserial=(float **)malloc((m+2*N-1)*sizeof(float *));
	for (i=0;i<(m+2*(N-1));i++)  {
		usingimageserial[i]=(float *)malloc((n+2*(N-1))*sizeof(float));
	}
	
	for (i=0;i<(m+2*(N-1));i++) {
		for (j=0;j<(n+2*(N-1));j++) {
			usingimageserial[i][j]=0;
		}
	}
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			usingimageserial[i+N-1][j+N-1]=imageserial[i][j];
		}
	}
	

	//kyrios vroxos kai xronometrhsh seiriakoy kwdika
	gettimeofday(&startwtime,NULL);
	for (i=0;i<level;i++)  {
		
		decompose(imageoutserial, filter, quadfilter, usingimageserial, i+1);
		
				
		for (k=0;k<((int)(m/pow(2.0,i+1))+2*(N-1));k++) {
			for (j=0;j<((int)(n/pow(2.0,i+1))+2*(N-1));j++) {
				usingimageserial[k][j]=0;
			}
		}
		for (k=0;k<(int)(m/pow(2.0,i+1));k++) {
			for (j=0;j<(int)(n/pow(2.0,i+1));j++) {
				usingimageserial[k+N-1][j+N-1]=imageoutserial[k][j];
			}
		}
		
		
		
	}
	gettimeofday(&endwtime,NULL);

	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);
	printf("Serial wall clock time = %f\n", seq_time);














	//memory allocation kai initialization twn pinakwn poy tha xrisimopoiithoun gia tin parallhlh ylopoiisi

	usingimage=(float *)malloc((m+2*(N-1))*(n+2*(N-1))*sizeof(float));
	
	for (i=0;i<m+2*(N-1);i++) {
		for (j=0;j<n+2*(N-1);j++) {
			usingimage[i*(n+2*(N-1))+j]=0;
		}
	}
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			usingimage[(i+N-1)*(n+2*(N-1))+j+N-1]=image[i*n+j];
		}
	}
	
	
	
	cudaMalloc(&devfilter,2*N*sizeof(float));
	cudaMemcpy(devfilter, filter, 2*N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc(&devquadfilter,2*N*sizeof(float));
	cudaMemcpy(devquadfilter, quadfilter, 2*N*sizeof(float), cudaMemcpyHostToDevice);
	griddim1=m*n/(2*a*b);
	
	sum1=0;
	sum2=0;
	sum3=0;
	for (i=0;i<level;i++) {
		int k1=(int) (m/pow(2.0,i));
		int k2=(int) (n/pow(2.0,i));
		sum1+=k1*(k2+2*(N-1));
		sum2+=(k1+2*(N-1))*k2;
	}
	int k1=(int) (m/pow(2.0,level));
	int k2=(int) (n/pow(2.0,level));
	sum1+=k1*(k2+2*(N-1));
	sum3=level*m*n;

	
	cudaMalloc(&devtempimage,sum1*sizeof(float));
	image1=(float *)malloc(sum1*sizeof(float));
	cudaMalloc(&devimageinterout,sum2*sizeof(float));
	image2=(float *)malloc(sum2*sizeof(float));
	cudaMalloc(&devimageout,sum3*sizeof(float));	
	image3=(float *)malloc(sum3*sizeof(float));

	
	for (i=0;i<sum1;i++) {
		image1[i]=0.0;
	}
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			image1[i*(n+2*(N-1))+j+N-1]=image[i*n+j];
		}
	}
	cudaMemcpy(devtempimage, image1, sum1*sizeof(float), cudaMemcpyHostToDevice);
	

	for (i=0;i<sum2;i++) {
		image2[i]=0.0;
	}
	cudaMemcpy(devimageinterout, image2, sum2*sizeof(float), cudaMemcpyHostToDevice);

	
	for (i=0;i<sum3;i++) {
		image3[i]=0.0;
	}
	cudaMemcpy(devimageout, image3, sum3*sizeof(float), cudaMemcpyHostToDevice);


	
	//kyrios vroxos kai xronometrhsh toy parallhloy kwdika
	gettimeofday(&startwtime,NULL);
	for (i=0;i<level;i++) {
		sum1=0;
		sum2=0;
		sum3=0;
		sum4=0;
		for (j=0;j<i;j++) {
			int k1=(int) (m/pow(2.0,j));
			int k2=(int) (n/pow(2.0,j));
			sum1+=k1*(k2+2*(N-1));
			sum2+=(k1+2*(N-1))*k2;
			sum3+=m*n;
		}
		int k1=(int) (m/pow(2.0,i));
		int k2=(int) (n/pow(2.0,i));
		sum4=sum1+k1*(k2+2*(N-1));
		dim3 grid1(griddim1, 2);
		dim3 block1(a, b);
		kernel1<<<grid1, block1>>>(devimageinterout+sum2, devtempimage+sum1, devfilter, devquadfilter, k1, k2, 2*N);


		dim3 grid2(griddim1/2, 4);
		dim3 block2(a, b);
		kernel2<<<grid2, block2>>>(devimageout+sum3, devtempimage+sum4, devimageinterout+sum2, devfilter, devquadfilter, k1, k2, n, 2*N);
		
		griddim1/=4;
	}
	gettimeofday(&endwtime,NULL);
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);
	printf("Parallel wall clock time = %f\n", seq_time);
	cudaMemcpy(image3, devimageout, level*m*n*sizeof(float), cudaMemcpyDeviceToHost);

	for (i=0;i<level;i++) {
		int k1=(int) (m/pow(2.0,i));
		int k2=(int) (n/pow(2.0,i));
		for (k=0;k<k1;k++) {
			for (j=0;j<k2;j++) {
				imageout[k*n+j]=image3[i*m*n+k*n+j];
			}
		}
	}

	//ypologismos megistou sfalmatos
	
	err=(float **)malloc(m*sizeof(float *));
	for (i=0;i<m;i++) {
		err[i]=(float *)malloc(n*sizeof(float));
	}
	
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			err[i][j]=imageoutserial[i][j]-imageout[i*n+j];
		}
	}
	
	
	
	

	float max;
	max=fabs(err[0][0]);
	
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			if (fabs(err[i][j])>max)
				max=fabs(err[i][j]);
		}
	}
	printf("\n   maximum absolute error = %f    \n", max);



	
}

//i synarthsh poy tha xrisimopoiithei gia thn aposynthesi tis eikonas ston seiriako kwdika
void decompose(float **target, float *filter, float *quadfilter, float **input, int level) {
	//memory allocation kai initialization twn pinakwn kai twn metavlhtwn poy tha xrisimopoiisei i decompose
	int end1=(int) (m/pow(2.0,level));
	int end2=(int) (n/pow(2.0,level));
	int i, j, k;
	float sum1, sum2, sum3, sum4;
	
	float **temptemp1, **temptemp2, **temp1, **temp2, **temp3, **temp4;
	temptemp1=(float **)malloc((m+2*(N-1))*sizeof(float *));
	for (i=0;i<m+2*(N-1);i++)  {
		temptemp1[i]=(float *)malloc((n+2*(N-1))*sizeof(float));
	}

	temptemp2=(float **)malloc((m+2*(N-1))*sizeof(float *));
	for (i=0;i<m+2*(N-1);i++)  {
		temptemp2[i]=(float *)malloc((n+2*(N-1))*sizeof(float));
	}


	temp1=(float **)malloc(m/2*sizeof(float *));
	for (i=0;i<m/2;i++)  {
		temp1[i]=(float *)malloc(n/2*sizeof(float));
	}
	
	temp2=(float **)malloc(m/2*sizeof(float *));
	for (i=0;i<m/2;i++)  {
		temp2[i]=(float *)malloc(n/2*sizeof(float));
	}

	temp3=(float **)malloc(m/2*sizeof(float *));
	for (i=0;i<m/2;i++)  {
		temp3[i]=(float *)malloc(n/2*sizeof(float));
	}

	temp4=(float **)malloc(m/2*sizeof(float *));
	for (i=0;i<m/2;i++)  {
		temp4[i]=(float *)malloc(n/2*sizeof(float));
	}
	
	
	//ypologismos toy endiamesou pinaka
	for (i=0;i<2*end1+2*(N-1);i++) {
		for (j=0;j<end2;j++) {
			sum1=0;
			sum2=0;
			for (k=0;k<2*N;k++) {
				
				sum1+=input[i][2*j+k]*filter[2*N-k-1];
				sum2+=input[i][2*j+k]*quadfilter[2*N-k-1];

			}
			temptemp1[i][j]=sum1;
			temptemp2[i][j]=sum2;
			
		}
	}
	
	//ypologismos toy telikoy pinaka
	for (j=0;j<end2;j++) {
		for (i=0;i<end1;i++) {
			sum1=0;
			sum2=0;
			sum3=0;
			sum4=0;
			
			for (k=0;k<2*N;k++) {
				
				sum1+=temptemp1[2*i+k][j]*filter[2*N-k-1];
				sum2+=temptemp1[2*i+k][j]*quadfilter[2*N-k-1];
				sum3+=temptemp2[2*i+k][j]*filter[2*N-k-1];
				sum4+=temptemp2[2*i+k][j]*quadfilter[2*N-k-1];
			}
			
			temp1[i][j]=sum1;
			temp2[i][j]=sum2;
			temp3[i][j]=sum3;
			temp4[i][j]=sum4;
		}
	}
	
	//eggrafi toy apotelesmatos
	for (i=0;i<end1;i++) {
		for (j=0;j<end2;j++) {
			target[i][j]=temp1[i][j];
			target[i][j+end2]=temp2[i][j];
			target[i+end1][j]=temp3[i][j];
			target[i+end1][j+end2]=temp4[i][j];
		}
	}
	
	free(temptemp1);
	free(temptemp2);
	free(temp1);
	free(temp2);
	free(temp3);
	free(temp4);
}


