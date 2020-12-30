*********************************
Tomer Locker
Digital Geometric Processing
*********************************

The code attached runs three diffrent parameterizations. From 3D to 2D using UV sets.

First you must attached the plugin to Maya. Go to Windows, Settings/Preferences, Plug-In Manager.
There you will see a DGP_CODE_DIR plugin. Check the Release Loaded box to attach the plugin.

In order to see the first two parameterizations, open an object of your choice with genus 0 (if its not it won't let you continue) 
and write the following command in the MEL window.

harmonicFlatteningCmd

Then in order to see the UV sets, right click on your object go down to UV Sets and you will see
Harmonic-Uniform
Harmonic-Cotangent

In order to see the third parameterization, open an object of your choice with genus 0 (if its not it won't let you continue) 
and write the following command in the MEL window.

LCSMCmd
Then in order to see the UV sets, right click on your object go down to UV Sets and you will see
LCSM

If at any point you dont see the UV set make sure your UV Editor is open by going to Pannels and UV Editor on the very bottom.

