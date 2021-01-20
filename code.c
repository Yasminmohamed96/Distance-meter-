
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "E:/programs/programming/keil/EE319Kware/inc/tm4c123gh6pm.h"
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

//Global variables 
int numbers[10] = {0x40,0x79,0x24,0x30,0x19,	// Each value turns on bits needed
		0x12,0x02,0x78,0x00,0x10}; // to show numbers in display
int digit1, digit2, digit3;		// Number to be displayed in each digit
uint32_t i =0;
uint32_t j;
void PortB_init(void)
{
    int delay; 
    SYSCTL_RCGC2_R |=2; //Port B
	  delay = 1; 
    GPIO_PORTB_LOCK_R= 0x4C4F434B;
    GPIO_PORTB_CR_R= 0xFF; 
    GPIO_PORTB_DIR_R=0x7F; // from PB0 to PB6 output 7seg,, and PB7 for input echo 
    GPIO_PORTB_DEN_R =0xFF;
    GPIO_PORTB_AFSEL_R =0x00;
    GPIO_PORTB_PCTL_R = 0x00000000;
    GPIO_PORTB_AMSEL_R = 0x00;
		GPIO_PORTB_PUR_R = 0x00;          // Enable pullup resistors on PB4,PF0  
}
void PortE_Init(void){
  uint32_t delay ; 
  SYSCTL_RCGC2_R |= 0x00000010;     // Port E clock initialized
  delay=1;
  GPIO_PORTE_LOCK_R = 0x4C4F434B;
  GPIO_PORTE_CR_R = 0x0F;           // Allow changes to PE3-0       
  GPIO_PORTE_AMSEL_R = 0x00;        // Disable analog function
  GPIO_PORTE_AMSEL_R = 0x00;
  GPIO_PORTE_PCTL_R = 0x00000000;   // GPIO clear bit PCTL  
  GPIO_PORTE_DIR_R = 0x0F;          // PE3-PE0 output          
  GPIO_PORTE_DEN_R = 0x0F;					// Enable digital pins PE3-PE0
	GPIO_PORTE_PUR_R = 0x00;          // Disable pullup resistors
}
void sys_tick_timer(void) //delay 10 micro
{
  NVIC_ST_CTRL_R=0;
  NVIC_ST_RELOAD_R=160; //(1*10)/0.0625
  NVIC_ST_CURRENT_R=0;
  NVIC_ST_CTRL_R =5;
  while((NVIC_ST_CTRL_R & 0x10000)==0){};
 }
  void SysTick_Init(void)
{
  NVIC_ST_CTRL_R=0;
  NVIC_ST_RELOAD_R=16; //(1)/0.0625 equivalent to 1 microsecond
  NVIC_ST_CURRENT_R=0;
  NVIC_ST_CTRL_R =5;
}
void SysTick_delay (uint32_t D)
{
	NVIC_ST_CTRL_R=0;
  NVIC_ST_CURRENT_R=0;
  NVIC_ST_CTRL_R =5;
	for(j=0;j<D;j=j+1)
{
  while((NVIC_ST_CTRL_R & 0x10000)==0){};
}
}
uint32_t Duration =0;
uint32_t Distance_CM =0;
 void systick_distance(void) // we can use it as a counter
{
  uint32_t ii=1;
  NVIC_ST_CTRL_R= 0;
  NVIC_ST_RELOAD_R = 16;
  NVIC_ST_CURRENT_R=0;
  NVIC_ST_CTRL_R =5;
  while((GPIO_PORTB_DATA_R&0x80)==0x80)
  {
	if((NVIC_ST_CTRL_R & 0x10000)==0x10000)
		ii++;
  }
	NVIC_ST_CTRL_R= 0;
  Duration = ((ii*NVIC_ST_RELOAD_R) - NVIC_ST_CURRENT_R)*0.0625; //in microsecond
  Distance_CM = (Duration*0.439/57.7);
	//if (Distance_CM<40)
//		Distance_CM=Distance_CM*0.448 - 2;	
	//else if (Distance_CM<100)
	//	Distance_CM=Distance_CM*0.448 - 4;
//	else if (Distance_CM>200&& Distance_CM<250)
	//	Distance_CM=Distance_CM*0.45 + 2 ;
//	else if (Distance_CM>250 && Distance_CM<300)
//		Distance_CM=Distance_CM*0.438 + 5;
//	else if (Distance_CM>300)
//		Distance_CM=Distance_CM*0.437 + 5;
//	else
	//	Distance_CM=Distance_CM*0.463 + 3;
 }  
 void systick_Distance(void) // we can use it as a counter
{
  NVIC_ST_CURRENT_R=0;
  while((GPIO_PORTB_DATA_R&0x80)==0x80)
  {
	while((NVIC_ST_CTRL_R & 0x10000)==0x00000)
		Duration++;
  }
  Distance_CM = (Duration/58);
 }
// Creates 0.1ms delay
void Delay2(void){
	unsigned long volatile time;
	time = 727240*200/91000;  // 0.1 ms
  while(time){
			time--;
  } 
}
// Takes digit and number for LED display
void Display(int digit, int number){
	GPIO_PORTB_DATA_R &= 0x80;		// Turns off LEDs
	GPIO_PORTE_DATA_R = digit;		// Selects digit
	GPIO_PORTB_DATA_R |= numbers[number];	// Turns on number in selected digit
	Delay2();					// Wait 0.1 ms
	 } 

// Splits number in counter into 4 separate numbers for each digit
void NumSplit(int counted){
	digit1 = counted%10;	//Copies value in counter, divides it by 10 and then keeps remainder
    counted /= 10;	//Dividing value in counter by 10 shifts it by one decimal
    digit2 = counted%10;
    counted /= 10;
    digit3 = counted%10;
}

int main (void){
   PortB_init ();
   PortE_Init();
   //SysTick_Init();
while (1){
				i=0;
	   GPIO_PORTE_DATA_R&=~0x08;//turn trigger off
				 sys_tick_timer ();
	   GPIO_PORTE_DATA_R|=0x08; //pin E3 IS trigger  
     	   sys_tick_timer (); //if flag reached one that mean 10 micro sec delay occured 
	   GPIO_PORTE_DATA_R&=~0x08;//turn trigger off
	   while ((GPIO_PORTB_DATA_R&0x80)==0){};//busy waiting until the echo pin is turned to 1 
		 systick_distance();
		while (i<1000)
		{		
		 NumSplit(Distance_CM);							// Split value in counter into 3 numbers
		 Display(1,digit1);								// Display number for lowest digit
		 Display(2,digit2);			
		 Display(4,digit3);
			 i=i+1;
		}
}
}
