/*
  pll.c

  This software accepts a bulk data transfer that contain register information
  for the ADF4351 PLL. This data is then output

 */
#include <stdio.h>

#include <fx2regs.h>
#include <fx2macros.h>
#include <serial.h>
#include <delay.h>
#include <autovector.h>
#include <lights.h>
#include <setupdat.h>
#include <eputils.h>


#define SYNCDELAY SYNCDELAY4
#define REARMVAL 0x80
#define REARM() EP2BCL=REARMVAL



volatile WORD bytes;
volatile __bit gotbuf;
volatile BYTE icount;
volatile __bit got_sud;
DWORD lcount;
__bit on;


//SPI interface definitions

#define SETDIN		IOA |= 0x04
#define CLRDIN		IOA &= 0xfb
#define SETCLK		IOA |= 0x02			// Pin assignments for TSM Proto Board
#define CLRCLK		IOA &= 0xfd
#define SETLE		IOA |= 0x01
#define CLRLE		IOA &= 0xfe
#define SETTRIG		IOA |= 0x09
#define CLRTRIG		IOA &= 0xf6

#define SETCE		IOA |= 0x40
#define CLRCE		IOA &= 0xBF

#define SETPDRF		IOA |= 0x20
#define CLRPDRF		IOA &= 0xDF


#define SETLDAC		IOA |= 0x10
#define CLRLDAC		IOA &= 0xef
#define SETCLR		IOA |= 0x08
#define CLRCLR		IOA &= 0xf7

// prototypes
void write_SPI(DWORD bits, DWORD d);
void set_ce_pdrf(int ce_pdrf_data);


void main() {

 DWORD r0;
 DWORD r1;
 DWORD r2;
 DWORD r3;
 DWORD r4;
 DWORD r5;

 REVCTL=0; // not using advanced endpoint controls

 d2off();
 on=0;
 lcount=0;
 got_sud=FALSE;
 icount=0;
 gotbuf=FALSE;
 bytes=0;

 // renumerate
 RENUMERATE_UNCOND();


 SETCPUFREQ(CLK_48M);
 SETIF48MHZ();
 sio0_init(57600);


 USE_USB_INTS();
 ENABLE_SUDAV();
 ENABLE_SOF();
 ENABLE_HISPEED();
 ENABLE_USBRESET();

 // IO port directions
 OEA = 0xff;						//direction bits on portA, All outputs
 IOA = 0x9f; 					//SYNC 1, SCLK 0, SDIN 0, LDAC 1, CLR 1, Bin2sC 1, RSTIN 1, CSADC 1
 OEB = 0x00;						//direction bits on portB, All inputs

 // IO port defaults
 CLRDIN;

 // only valid endpoints are 2/6
 EP2CFG = 0xA2; // 10100010
 SYNCDELAY;
 EP6CFG = 0xE2; // 11100010
 SYNCDELAY;
 EP1INCFG &= ~bmVALID;
 SYNCDELAY;
 EP1OUTCFG &= ~bmVALID;
 SYNCDELAY;
 EP4CFG &= ~bmVALID;
 SYNCDELAY;
 EP8CFG &= ~bmVALID;
 SYNCDELAY;

 // arm ep2
 EP2BCL = 0x80; // write once
 SYNCDELAY;
 EP2BCL = 0x80; // do it again


 EA=1; // global interrupt enable
 printf ( "Done initializing stuff\n" );

 // enable PLL chip and turn off RF mute
 set_ce_pdrf(0x1C);

 d3off();

 // activate default register settings

// set one : 1 GHz PLL
//r0=0x00C80000;
//r1=0x08008011;
//r2=0x00004E42;
//r3=0x000004B3;
//r4=0x00A5003C;
//r5=0x00580005;

 // set two : 125 MHz PLL
 r0 = 0x00C80000;
 r1 = 0x08008011;
 r2 = 0x00004E42;
 r3 = 0x000004B3;
 r4 = 0x00D5003C;
 r5 = 0x00580005;

 write_SPI(32,r5);
 delay(10);
 write_SPI(32,r4);
 delay(10);
 write_SPI(32,r3);
 delay(10);
 write_SPI(32,r2);
 delay(10);
 write_SPI(32,r1);
 delay(10);
 write_SPI(32,r0);
 delay(10);

 while(TRUE) {

  if ( got_sud ) {
      printf ( "Handle setupdata\n" );
      handle_setupdata();
      got_sud=FALSE;
  }

// Data is received in fifo 2 and returned in fifo 6

  if ( !(EP2468STAT & bmEP2EMPTY) ) {
       printf ( "ep2 out received data\n" );
      if  ( !(EP2468STAT & bmEP6FULL) ) { // wait for at least one empty in buffer
                 WORD i;
                 printf ( "Sending data to ep6 in\n");

                 bytes = MAKEWORD(EP2BCH,EP2BCL);

                 // fifo data holds the PLL registers
                 // order is r5,r4,r3,r2,r1,r0 (4 bytes per register)
                 // MSB to LSB
                 r5 = EP2FIFOBUF[3] + EP2FIFOBUF[2]*256 + EP2FIFOBUF[1]*65536 + EP2FIFOBUF[0]*16777216;
                 r4 = EP2FIFOBUF[7] + EP2FIFOBUF[6]*256 + EP2FIFOBUF[5]*65536 + EP2FIFOBUF[4]*16777216;
                 r3 = EP2FIFOBUF[11] + EP2FIFOBUF[10]*256 + EP2FIFOBUF[9]*65536 + EP2FIFOBUF[8]*16777216;
                 r2 = EP2FIFOBUF[15] + EP2FIFOBUF[14]*256 + EP2FIFOBUF[13]*65536 + EP2FIFOBUF[12]*16777216;
                 r1 = EP2FIFOBUF[19] + EP2FIFOBUF[18]*256 + EP2FIFOBUF[17]*65536 + EP2FIFOBUF[16]*16777216;
                 r0 = EP2FIFOBUF[23] + EP2FIFOBUF[22]*256 + EP2FIFOBUF[21]*65536 + EP2FIFOBUF[20]*16777216;

                 write_SPI(32,r5);
                 delay(10);
                 write_SPI(32,r4);
                 delay(10);
                 write_SPI(32,r3);
                 delay(10);
                 write_SPI(32,r2);
                 delay(10);
                 write_SPI(32,r1);
                 delay(10);
                 write_SPI(32,r0);
                 delay(10);

                 // loopback to verify transfer

                 for (i=0;i<bytes;++i) EP6FIFOBUF[i] = EP2FIFOBUF[i];

                 // can copy whole string w/ autoptr instead.
                 // or copy directly from one buf to another

                 // ARM ep6 out
                 EP6BCH=MSB(bytes);
                 SYNCDELAY;
                 EP6BCL=LSB(bytes);

                 REARM(); // ep2
                 //printf ( "Re-Armed ep2\n" );

         }
   }
 }

}

// copied routines from setupdat.h

BOOL handle_get_descriptor() {
  return FALSE;
}

// value (low byte) = ep
#define VC_EPSTAT 0xB1

BOOL handle_vendorcommand(BYTE cmd) {

 switch ( cmd ) {

     case VC_EPSTAT:
        {
         __xdata BYTE* pep= ep_addr(SETUPDAT[2]);
         printf ( "ep %02x\n" , *pep );
         if (pep) {
          EP0BUF[0] = *pep;
          EP0BCH=0;
          EP0BCL=1;
          return TRUE;
         }
        }
     default:
          printf ( "Need to implement vendor command: %02x\n", cmd );
 }
 return FALSE;
}

// this firmware only supports 0,0
BOOL handle_get_interface(BYTE ifc, BYTE* alt_ifc) {
 printf ( "Get Interface\n" );
 if (ifc==0) {*alt_ifc=0; return TRUE;} else { return FALSE;}
}
BOOL handle_set_interface(BYTE ifc, BYTE alt_ifc) {
 printf ( "Set interface %d to alt: %d\n" , ifc, alt_ifc );

 if (ifc==0&&alt_ifc==0) {
    // SEE TRM 2.3.7
    // reset toggles
    RESETTOGGLE(0x02);
    RESETTOGGLE(0x86);
    // restore endpoints to default condition
    RESETFIFO(0x02);
    EP2BCL=0x80;
    SYNCDELAY;
    EP2BCL=0X80;
    SYNCDELAY;
    RESETFIFO(0x86);
    return TRUE;
 } else
    return FALSE;
}

// get/set configuration
BYTE handle_get_configuration() {
 return 1;
 }
BOOL handle_set_configuration(BYTE cfg) {
 return cfg==1 ? TRUE : FALSE; // we only handle cfg 1
}


// copied usb jt routines from usbjt.h
void sudav_isr() __interrupt SUDAV_ISR {

  got_sud=TRUE;
  CLEAR_SUDAV();
}

__bit on5;
__xdata WORD sofct=0;
void sof_isr () __interrupt SOF_ISR __using 1 {
    ++sofct;
    if(sofct==8000) { // about 8000 sof interrupts per second at high speed
        on5=!on5;
        if (on5) {d5on();} else {d5off();}
        sofct=0;
    }
    CLEAR_SOF();
}

void usbreset_isr() __interrupt USBRESET_ISR {
    handle_hispeed(FALSE);
    CLEAR_USBRESET();
}
void hispeed_isr() __interrupt HISPEED_ISR {
    handle_hispeed(TRUE);
    CLEAR_HISPEED();
}


/* PLL chip and RF enable combo */
void set_ce_pdrf(int ce_pdrf_data)
{
      CLRCE;
	  CLRPDRF;

 	if(ce_pdrf_data == 0x4)
    {
		CLRCE;
		CLRPDRF;
	}

	if(ce_pdrf_data == 0x14)
    {
		CLRCE;
		SETPDRF;
	}

	if(ce_pdrf_data == 0xC)
    {
		SETCE;
		CLRPDRF;
	}

    if(ce_pdrf_data == 0x1C)
    {
		SETCE;
		SETPDRF;
	}

}

/* PLL write serial word using LE based control */
void write_SPI(DWORD bits, DWORD d)
{
DWORD mask;
DWORD j;
BYTE bit;

	//BitTest = (pow (2,Bits)) / 2;
  mask = 1 << (bits -1);

  CLRDIN;
	CLRLE;

	for(j = 0; j < bits; j++)	{

      bit = !((d & mask) == 0);

      if (bit) {
				SETDIN;
      }	else {
				CLRDIN;
      }

      mask = mask >> 1;

      delay(1);
			SETCLK;
      delay(1);
			CLRCLK;
      delay(1);
  }
  delay(1);
	SETLE; // latch data to active
  delay(1);
	CLRLE;
  CLRDIN;

}
