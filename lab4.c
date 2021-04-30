#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stm32l4xx_hal.h"
#include "stm32l476g_discovery.h"
#include "arm_math.h"
#include "ece486.h"
#define PI 3.14159265

extern FlagStatus KeyPressed;
int main()
{
	int num_of_samps = 2048; //need to figure out this value
	int Fo = 0; // hz value
	double f0 = 0; //cycles per sample
	double angle; //variable to hold the angle for the cosine and sin calls
	char lcd_str[8];
	
	//set the number of coefficents and the fir coefficients from the output of the matlab filter designer. 
	int n_coefs = 13;
	float fir_coefs[13] = {
  	-0.004474141829264,  0.01741756658841, -0.04283986607915,  0.07971359883835,
  	  -0.1204124427275,   0.1526398872318,   0.8350608095554,   0.1526398872318,
    	-0.1204124427275,  0.07971359883835, -0.04283986607915,  0.01741756658841,
  	-0.004474141829264
	};
	//declare the two instances needed for the decimate function and cfft function. 
	arm_fir_decimate_instance_f32 *filter;
	arm_cfft_instance_f32 *cfft;

	initialize_ece486(FS_48K, MONO_IN, STEREO_OUT, MSI_INTERNAL_RC); //set up all parts of stm board from hummels example code

	//declare all arrays needed for testing
	float input, sweeping_freq, mag_calced, real_part, imag_part, real_decimate, imag_decimate, combined_decimate;
	
	//alocate memory for all arrays used based on number of samples
	input = (float*)calloc(num_of_samps * 10,sizeof(float));
	sweeping_freq = (float*)calloc(num_of_samps,sizeof(float));
	mag_calced = (float*)calloc(num_of_samps,sizeof(float));
	real_part = (float*)calloc(num_of_samps * 10,sizeof(float));
	imag_part = (float*)calloc(num_of_samps * 10,sizeof(float));
	real_decimate = (float*)calloc(num_of_samps,sizeof(float));
	imag_decimate = (float*)calloc(num_of_samps,sizeof(float));
	combined_decimate = (float*)calloc(num_of_samps*2,sizeof(float));// 2 times the size as rest to hold two arrays
	
	//declare state and allocate memory for state
	float *state = (float *)malloc((n_coef+num_of_samps-1)*sizeof(float));//allocate memory for based on blocksize and number of coefficients
	if (state == NULL) {
		flagerror(MEMORY_ALLOCATION_ERROR);
		while(1);
	}

	//check to make sure all memory was allocated correctly
	if (input == NULL || sweeping_freq == NULL 
	|| mag_calced == NULL || real_part == NULL 
	|| imag_part == NULL || real_decimate == NULL 
	|| imag_decimate == NULL || combined_decimate  || state== NULL ||) {
		flagerror(MEMORY_ALLOCATION_ERROR);
		while(1);
	}
	for (int k = 0; k < num_of_samps; k++) {
		sweeping_freq[i] = -1 + (i/1024); //set the sweeping frequncy to a range from -1 to 1
	}

	arm_fir_decimate_init_f32(filter,n_coefs,25,fir_coeffs,state,num_of_samps); //initialize the fir filter used in the decimation function 
	arm_cfft__init_f32(cfft,num_of_samps); //initialize the cfft instance used in the cfft function

	while(1){ //forever loop
		
		getblock(input);
		f0 = F0 / 32; //update fo value based on 
		for(int j = 0; j < num_of_samps * 10; j++) {
			angle = (2*PI*f0*j) % (2 * PI)
			real_part[i] = input[i] * arm_cos_f32(angle); find the real part of the input
			imag_part[i] = input[i] * arm_sin_f32(angle) * -1; //find the imaginary part of the input
		}
		//find real filter value and thats mostly it
		arm_fir_decimate_f32(filter, real_part, real_decimate, num_of_samps); //use decimate function to filter and ecimate the real portion of the input
		arm_fir_decimate_f32(filter, imag_part, imag_decimate, num_of_samps); //use decimate function to filter and decimate the imaginary portion of the output
		
		// combine the real and imaginary decimation portions
		for(int i = 0; i < num_of_samps * 2; i+2) {
			combined_decimate[i] = real_decimate[i];
			combined_decimate[i+1] = imag_decimate[i];
		}
		
		// use the fast fourier transform on the combined decimation array. 
		arm_cfft_f32(cfft,combined_decimate,0,0);
		
		// find the complex magnitude at every output place these values in out_2
		arm_cmplex_mag_f32(combined_decimate,mag_calced,num_of_samps);
		
		//scale the output values to meet specification
		for (int l = 0; l < num_of_samps; l++) {
			mag_calc[i] = mag_calced[i] * .011718 * 3; //scale the magnitude by the resoultion and the max output of the DAC. 
		}

		// output 1 and output 2 to DAC with hummels example
		putblockstereo(sweeping_freq,mag_calced);

		if (KeyPressed) { // check if key is pressed
			KeyPressed = RESET; //reset the key
			if (F0 = 15) { //if F0 is already at max, put back to min
				F0 = 0;
			}
			else { //else add one to F0
				F0++;
			}
		}
		sprintf(lcd_str, "F0: %2dkHz",F0); // place F0 value in the lcd string
		BSP_LCD_GLASS_DisplayString( (uint8_t *)lcd_str); //print the F0 value to the display
	}

	return 0; //end of program
}}
