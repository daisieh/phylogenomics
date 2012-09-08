RequireVersion ("0.9920060830");
VERBOSITY_LEVEL = -1;

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

ModelMatrixDimension = 64;
for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k] == 10)
	{
		ModelMatrixDimension = ModelMatrixDimension -1;
	}
}

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"2RatesAnalyses"+DIRECTORY_SEPARATOR+"MG94xREV.mdl");

SetDialogPrompt     ("Choose a nucleotide alignment");
DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);

DataSetFilter	  	filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

SKIP_MODEL_PARAMETER_LIST = 1;
done 					  = 0;

while (!done)
{
	fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
	fscanf  (stdin,"String", modelDesc);
	if (Abs(modelDesc)==6)
	{
		done = 1;
	}
}
modelType 				  = 0;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MG94custom.mdl");
SKIP_MODEL_PARAMETER_LIST = 0;

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");

brNames				= BranchName (givenTree,-1);
COVARIANCE_PARAMETER 				= {};
global global_OMEGA = 1;
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	ExecuteCommands ("givenTree."+brNames[k]+".nonSynRate:=givenTree."+brNames[k]+".omega*givenTree."+brNames[k]+".synRate;");
	COVARIANCE_PARAMETER["givenTree."+brNames[k]+".omega"] = 1;
}

LikelihoodFunction  theLnLik = (filteredData, givenTree);


for (k=0; k < Columns (brNames)-1; k=k+1)
{
//  set all of the branches to have omega constrained to the same global_OMEGA
	ExecuteCommands ("givenTree."+brNames[k]+".omega:=global_OMEGA;");
}

fprintf 					   (stdout, "\nFitting the global model to the data...\n");
Optimize 					   (res_global, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

for (k=0; k < Columns (brNames)-1; k=k+1)
{
//  set each branch to have unconstrained omega
	ExecuteCommands ("givenTree."+brNames[k]+".omega=global_OMEGA;");
}

fprintf 					   (stdout, "\nFitting the local model to the data...\n");
Optimize 					   (res_local, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

for (k=0; k < Columns (brNames)-1; k=k+1)
{
//  set each branch to have omega = 1
	ExecuteCommands ("givenTree."+brNames[k]+".omega:=1;");
}

fprintf 					   (stdout, "\nFitting the neutral model to the data...\n");
Optimize 					   (res_neutral, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

LR = 2(res_local[1][0]-res_global[1][0]);
DF = res_local[1][1]-res_global[1][1];

fprintf (stdout, "\nGlobal omega calculated to be ", R, "\n");

fprintf (stdout, "\nLRT for variable omega across the tree\n\tLR = ", LR, "\n\tConstraints = ", DF, "\n\tp-value = ", 1-CChi2(LR,DF),"\n\n");

LR = 2(res_global[1][0]-res_neutral[1][0]);
DF = res_global[1][1]-res_neutral[1][1];

fprintf (stdout, "\nLRT for single omega across the tree\n\tLR = ", LR, "\n\tConstraints = ", DF, "\n\tp-value = ", 1-CChi2(LR,DF),"\n\n");

COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (covMx, theLnLik);
//
VERBOSITY_LEVEL = 0;
//
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	fprintf (stdout, "Branch :", brNames[k], "\n\tomega MLE = ", Format (covMx[k][1],6,3), "\n\t95% CI = (",Format (covMx[k][0],6,3), ",", Format (covMx[k][2],6,3), ")\n");
}
