/*Visualization #1: Love plot*/

/*Naming the library where the trial data are being kept*/

LIBNAME cohort "C:\Users\webst\OneDrive\Documents\SER seminar data\Trial";

/*Naming the library where the target data are being kept*/

LIBNAME targ "C:\Users\webst\OneDrive\Documents\SER seminar data\Target";

/*First, let's bring the trial data into the "work" directory*/

DATA trialset;
	SET cohort.analytictrialcohort;
RUN;

/*Repeating something similar for the target data*/

DATA targset;
	SET targ.analytictargcohort;
RUN;

ods listing gpath= "C:\Users\webst\OneDrive\Documents\SER seminar data\Visualizations"; /*Assigning folder where figures will be saved*/

%macro Loveplot(TRIALDATA=, /*Name of trial/study population file*/
				TARGETDATA=, /*Name of target population file*/
				FIGURENAME=/*Specify the figure name*/);

/*First, we should give the trial and target variables indicating their origin*/

DATA _trial1;
	SET &TRIALDATA;
	intrial = 1; /*Giving the study/trial participants a variable indicating their origin*/
RUN;

DATA _targ1;
	SET &TARGETDATA;
	intrial = 0; /*Giving the target participants a variable indicating their origin*/
RUN;

/*Next, we concatenate the study/trial and target data*/

DATA _concat;
	SET _trial1 _targ1;
RUN;

/*And we create data sets including the materials we need for the SMD: for continuous variables,
we need the mean in the trial, the mean in the target, and the overall standard deviation, while for categorical ones,
we need the proportion in the trial and the proportion in the target.*/

/*Let's start with the continuous variable: age*/

/*Starting with obtaining the trial mean*/

PROC MEANS DATA=_concat NOPRINT;
	VAR age;
	WHERE intrial=1; /*Restricting to trial members*/
	OUTPUT OUT=_trialmean mean=mean_trial stddev=stddev_trial /*Creating an output dataset*/;
RUN;

/*Now the target mean*/

PROC MEANS DATA=_concat NOPRINT;
	VAR age;
	WHERE intrial=0; /*Restricting to target members*/
	OUTPUT OUT=_targmean mean=mean_targ stddev=stddev_targ /*Creating an output dataset*/;
RUN;

DATA _trialmeanage;
	LENGTH varname $ 20.;
	SET _trialmean;
	varname = "age";
RUN;

DATA _targmeanage;
	LENGTH varname $ 20.;
	SET _targmean;
	varname = "age";
RUN;


/*Now we combine those data sets*/

DATA _mergedage;
	MERGE _trialmeanage _targmeanage;
	BY varname;
RUN;
	
/*Next, we can create a macro to automate this task for the binary variables*/

%macro prepcatsmd(varname=,);

/*Creating a data set with the prevalence in the trial*/

PROC FREQ DATA=_concat NOPRINT;
	WHERE intrial=1;
	TABLES &varname/ OUT=_trialprev&varname /*creating the new data set*/;
RUN;

/*Creating a data set with the prevalence in the target*/

PROC FREQ DATA=_concat NOPRINT;
	WHERE intrial=0;
	TABLES &varname/ OUT=_targprev&varname /*creating the new data set*/;
RUN;

/*Altering these datasets to format them properly for joining*/

DATA _trialprev&varname;
	SET _trialprev&varname;
	LENGTH varname $ 20.;
	WHERE &varname = 1; /*Limiting to the row with covariate values of 1*/
	trialprev = PERCENT/100; /*Creating a new variable trialprev based on the default PERCENT output by SAS*/
	varname = "&varname"; /*Creating a variable varname with the same formatting as the one we did in the age setup and assigning it the correct value*/
	KEEP varname trialprev;
RUN;

DATA _targprev&varname;
	SET _targprev&varname;
	LENGTH varname $ 20.;
	WHERE &varname = 1; /*Limiting to the row with covariate values of 1*/
	targprev = PERCENT/100; /*Creating a new variable trialprev based on the default PERCENT output by SAS*/
	varname = "&varname"; /*Creating a variable varname with the same formatting as the one we did in the age setup and assigning it the correct value*/
	KEEP varname targprev;
RUN;

DATA _merged&varname;
	MERGE _trialprev&varname _targprev&varname;
	BY varname;
RUN;

%MEND;

/*Calling the macro for each binary variable*/

%prepcatsmd(varname=female);
%prepcatsmd(varname=livermets);
%prepcatsmd(varname=colon);
%prepcatsmd(varname=WildKRAS);

/*Combining all these data sets into one dataset with what we need to calculate all the covariate SMDs*/

DATA SMDs;
	SET _mergedage _mergedlivermets _mergedWildKRAS _mergedfemale _mergedcolon;
	IF varname = "age" THEN DO; /*Calculating SMDS for age*/
		SMD = ( mean_trial - mean_targ ) / ( ( ( stddev_trial ** 2 + stddev_targ ** 2) / 2 ) ** 0.5 );
	END;
	ELSE DO; /*Calculating SMDS for binary variables*/
		SMD = ( trialprev - targprev ) / (  ( ( trialprev * (1 - trialprev) + targprev * (1 - targprev) ) / 2 ) ** 0.5 );
	END;
RUN;



ODS GRAPHICS / reset=index imagename="&FIGURENAME"; /*Assigning the image name*/

/*Creating the figure*/

/*Cleaning up the figure*/

DATA ForplotSMDs;
	SET SMDs;
	LENGTH varname2 $ 25.;
	IF varname="colon" THEN varname2 = "Colon cancer";
	IF varname="female" THEN varname2 = "Female sex";
	IF varname="livermets" THEN varname2 = "Liver metastases";
	IF varname="WildKRAS" THEN varname2 = "Wild-type KRAS";
	IF varname="age" THEN varname2 = "Age (in years)";
RUN;

PROC SGPLOT DATA=ForPlotSMDs;
	SCATTER Y=varname2 X=SMD / markerattrs=(color=black symbol=circlefilled);
	REFLINE 0 / AXIS=X lineattrs=(color = black pattern=solid);
	REFLINE 0.1 -0.1 / AXIS=X lineattrs=(color = gray pattern=dash);
	YAXIS LABEL="Covariate" discreteorder=data LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABEL="Standardized mean difference" MIN=-0.7 MAX=0.7 LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;

%mend;

/*Executing the macro*/;

%Loveplot(TRIALDATA=trialset,
		TARGETDATA=targset,
		FIGURENAME=Love_plot);
