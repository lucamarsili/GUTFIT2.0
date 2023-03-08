# gutfit

GUTFIT is a program which scan the parameter space of a GUT model in order to fit fermion masses and mixing parameters. 

The gutfit folder contains, typeIandIIseewsaw_v4.py which contains the info about the model used. It contains also experimentalneutrinomassmatrix.py which contains the data parameters used for the fit. 

Enter in the main folder $GUTFIT_-gutfit_chisquared/
To run the complete scan run: 
 python3 examples/multinest_v4.py examples/param_card3sigma_v4.dat -o $OUTPUT_FOLDER_NAME
 
To compute the measure of a certain point in the parameter space run: 
python3 examples/multinest_v4.py examples/param_card3sigma_v4.dat -o $OUTPUT_FOLDER_NAME $NAME_OF_PARAMETER_FILE.txt $NAMEPLOT.pdf

An example for the parameter file of the chosen benchmark point is BP1_220900021.txt


