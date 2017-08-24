/*
  modified from State change detection (edge detection)
 	
 Often, you don't need to know the state of a digital input all the time,
 but you just need to know when the input changes from one state to another.
 For example, you want to know when a button goes from OFF to ON.  This is called
 state change detection, or edge detection.
 
 This example shows how to detect when a button or button changes from off to on
 and on to off.
 	
 The circuit:
 * pushbutton attached to pin 2 from +5V
 * 10K resistor attached to pin 2 from ground
 * LED attached from pin 13 to ground (or use the built-in LED on
 most Arduino boards)
 
 The Arduino board contains a 6 channel (8 channels on the Mini and Nano, 16 on the Mega), 
 10-bit analog to digital converter. This means that it will map input voltages
 between 0 and 5 volts into integer values between 0 and 1023. This yields a 
 resolution between readings of: 5 volts / 1024 units or, .0049 volts (4.9 mV) 
 per unit. The input range and resolution can be changed using analogReference().
 
 */

// this constant won't change:
const int  SCTriggerAnalogInPin = A0;    // the pin that the laser signal from soundcard is attached to
const int  LaserAnalogInPin = A1;    // the pin that the laser signal from soundcard is attached to
const int SCTriggerOutPin = 0;       // pin to drive laser
const int LaserOutPin = 1;       // pin to drive laser
const int threshold = 10;       // 10 AD units is ~50 mV

// Variables will change:
int SCTriggerValue = 0;  // variable to store the value coming from the sensor
int LaserValue = 0;  // variable to store the value coming from the sensor

void setup() {
  // initialize the LED as an output:
  pinMode(SCTriggerOutPin, OUTPUT);
  pinMode(LaserOutPin, OUTPUT);
}


void loop() {
  // read the value from the sensor:
  SCTriggerValue = analogRead(SCTriggerAnalogInPin);
  if (SCTriggerValue>threshold) {
    digitalWrite(SCTriggerOutPin, HIGH);
  }
  else if (SCTriggerValue<=threshold) {
    digitalWrite(SCTriggerOutPin, LOW);
  }

  LaserValue = analogRead(LaserAnalogInPin);
  if (LaserValue>threshold) {
    digitalWrite(LaserOutPin, HIGH);
  }
  else if (LaserValue<=threshold) {
    digitalWrite(LaserOutPin, LOW);
  }






}













