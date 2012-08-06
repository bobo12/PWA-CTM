!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   THIS CD-ROM contains all the information of the EEC291 final PROJECT
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Authors : Boris Prodhomme
	  Pierre-Antoine Benard

Project Title : Implementation of a density based cell transmission model for
traffic estimation on highway using Ensemble Kalman Filter

!----------------------------------------------------------------------------!

The main matlab script is mainScript.m
In order to be able to run it, one needs an access to the Mobile Millennium database.
All the scripts must be run with a recent version of Matlab

Organization of the CDROM :

	- Literature folder : contains all the reading we used for the theory 
			      and implementation.

	- Final_Presentation : contains the Prezi file used for the presentation.
		To open it : 
			on windows : simply press the Prezi.exe app, it will open
				 a flash window in which will be able to navigate
			     	 through the presentation
			on Macs : similarly just launch the prezi file. (compatibility OK)
	
	- Simulation_Results : contains .png and .fig files from matlab simulations 
				(on the validation of both the CTM and the EnKF). 
				Some of the figures are used in the report.
	
	- Report : Contains all the source files necessary to generate the .pdf report from 
			the .tex. The subfolder figures contains mainly figures from the simulation
			on real data  and some extra flow charts for the report.
	- codes : 
		- data: csv file containing measurements data on Highway-I880 2012/03/05 Monday
		- sql: sql queries obtaining data
		- matlab: all the scripts
			Test_Enkf.m compares analytical solution and output of the Enkf
			- codebeta (all the codes that generate the figures verifying the validity
			of the CTM script )
				Example_Verify_CTM.m
				Verify_CTM_Daganzo.m are the two scripts that can be run
			- OtherEnkf (other implementations of the Enkf algorithms, they can be more 
			efficient in some cases)
			
			
			SCRIPTS TO RUN TO GET THE FIGURES
			-	Without access to the database :
				TestRouteI880.m (with Greenshields fundamental diagram)
				TestRouteI880_Daganzo.m
		
			- With access to the Mobile Millennium database
				mainScript.m (This script is completely automated, you just have 
				to insert the parameters and everything is done for you)


Enjoy!
