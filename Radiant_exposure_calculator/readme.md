rad_exposure_calc.m calculates the radiant exposure of an incident pulsed laser beam on your sample in the commonly used unit Milli-Joules per square centimeters.
This .m script was created in the context of pump probe spectroscopy research, so you'll find variables for both pump and probe radiant exposure. They are of course calculated in the same way.

You must provide 

* angle of incidence (in °, degree)
* power of your beam right at the spot on the sample (in mW)
* the size of your pump and probe spots on the sample. Povide it with left and right boundaries in the two orthogonal dimensions x and y of your beam profile. Keep in mind: In the resulting calculation, the beam shape in the surface plane is assumed as perfectly round, so you might have a deviation when your shape is oval.
* repetition rate of your laser light (in Hz)
