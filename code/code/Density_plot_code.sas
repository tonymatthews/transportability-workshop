/*Visualization #3: A density plots of the predicted "sampling" probabilities and their overlap*/

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


/*Next, getting the quantities we need for the chart*/;

%macro Densplot(	TRIALDATA=, /*Name of trial/study population file*/
				TARGETDATA=, /*Name of target population file*/
				CATPARAMS=,/*The categorical variables for the model*/
				MODELPARAMS=,/*All the variables to be used in the model*/
				FIGURENAME=);

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
	BY replicate;
RUN;

/*And we run the same logistic regression we used to to get the sampling ORs; this time, instead
of outputting the model parameters, we output the predicted probabilities to a data set named "withprobs"*/

PROC GENMOD DATA=_concat DESCENDING;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS / link=logit dist=binomial;
	OUTPUT OUT=_withprobs prob=sampprob;
RUN;


/*To make graphing easy, we can create two new variables-one for those in the trial, one for those
in the target*/

DATA _fordensityplot;
	SET _withprobs;
	IF intrial = 1 THEN trialprobs = sampprob;
	IF intrial = 0 THEN targprobs = sampprob;
RUN;

ODS GRAPHICS / reset=index imagename="&FIGURENAME"; /*Assigning the image name*/

/*Now, we plot histograms*/

PROC SGPLOT DATA=_fordensityplot NOAUTOLEGEND;
	HISTOGRAM trialprobs / fillattrs=(color=blue) transparency = 0.5 BINDWIDTH = 0.005;
	HISTOGRAM targprobs /  fillattrs=(color=red) transparency = 0.5 BINDWIDTH = 0.005;
	XAXIS LABEL="Predicted probability of presence in the trial data set" MIN=0 LABELATTRS=(size=14) VALUEATTRS=(size=12);
	YAXIS LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;

%mend;

/*Executing the macro*/;

%Densplot(TRIALDATA=trialset,
		TARGETDATA=targset,
		CATPARAMS=female colon livermets WildKRAS,
		MODELPARAMS=age female colon livermets WildKRAS,
		FIGURENAME=Density_plot);
