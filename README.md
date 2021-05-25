# sharc_cobramm_interface
This repository contains the files to run QM/MM trajectory surface hopping interfacing the SHARC and COBRAMM software. 

The are three kinds of files:

1) auxiliary scripts: i) setup_cobramm_single_point.py ii) setup_cobramm_init.py, iii) setup_cobramm_traj.py iv) combine_initconds.py

2) SHARC interfaces: updated i) SHARC_RICC2.py ii) SHARC_ORCA.py iii) SHARC_MOLCAS.py that include the option for cobramm and the ad hoc interface iv) SHARC_COBRAMM.py

3) sharc_cobramm_documentation.pdf where information about the interface, the two software and the other files listed before is reported.

Additional requirements are: both python 2.x and 3.y installed, SHARC (version 2.1) and COBRAMM (version 2.3), AMBER (version tested 17,18,19,20), a QM code among TURBOMOLE, ORCA, MOLCAS, OpenMOLCAS installed.

For any further information, consult the enclosed documentation or contact the development teams for questions, feedback or bug reports.   

