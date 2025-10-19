Theory: QCQP to SOCP

1. Proof using mostly algebra 
2. Use epigraph transform as mentioned in slides. 
   Make constraint SOCP compliant through rotated second order cone (was a bit more challenging, got some hints from GenAI)

Done by: Benjamin Labrecque

Theory: Linear Programming -- transform

Use epigraph transform again. Use other transformation from slides to convert constraints to requested form. 

Done by: Benjamin Labrecque

Theory: Linear Programming -- general LP to standard form

Done by: Horacio Lisdero

Programming Exercise: Mass Spring System

1. Solution Description

- f, grad, hessian: compute by hand and implement as code. Some hints from GenAI to make 
hessian more compact. 
- graph setup: add springs to neighbours if they are within bounds
- assemble gloabl energy, grad, hessian: use exercise description and video on course website
to understand how to assemble the global values from the local ones. 

2. Encountered Problems and Solutions

- Some tests fail when using precise equality: update them to use approx equality. 
- Some tests were a bit confusing (see Ilias discussion)

3. Done by: Benjamin Labrecque
