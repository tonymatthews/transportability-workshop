/*Visualization #4: A skyscraper plot of the weights resulting from the sampling model*/


/*Naming the library where the trial data are being kept*/

LIBNAME cohort "/local/users/mawc/PCORI/SER Seminar/Data/Trial data";

/*Naming the library where the target data are being kept*/

LIBNAME targ "/local/users/mawc/PCORI/SER Seminar/Data/Target data";


/*First, let's bring the trial data into the "work" directory*/

DATA trialset;
	SET cohort.analytictrialcohort;
RUN;

/*Repeating something similar for the target data*/

DATA targset;
	SET targ.analytictargcohort;
RUN;

ods listing gpath= "/local/users/mawc/PCORI/SER Seminar/Visualizations"; /*Assigning folder where figures will be saved*/

/*Next, getting the quantities we need for the chart*/;

%macro Skyscraper(	TRIALDATA=, /*Name of trial/study population file*/
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
RUN;

/*And we run the same logistic regression we used to to get the sampling ORs; this time, instead
of outputting the model parameters, we output the predicted probabilities to a data set named "withprobs"*/

PROC GENMOD DATA=_concat DESCENDING;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS / link=logit dist=binomial;
	OUTPUT OUT=_withprobs prob=sampprob;
RUN;

/*This time, we also fit an empty logit model to help us stabilize the N after weighting (meaning
the trial will have its original number of participants)*/

PROC GENMOD DATA=_withprobs DESCENDING;
	MODEL intrial =  / link=logit dist=binomial;
	OUTPUT OUT=_forweights prob=numerator;
RUN;

/*Now, we use this "forweights" data set to calculate stabilized weights for each participant*/

DATA _withweights;
	SET _forweights;
	IF intrial = 0 THEN DO; /*First setting every non-participant's weight to 0*/
		weight = 0;
	END;
	ELSE IF intrial = 1 THEN DO;
		weight = (numerator / (1 - numerator) ) / (sampprob / (1 - sampprob ) ); /*Calculating the stabilized inverse odds weight by transforming the probabilities*/
	END;
RUN;

/*Creating a new dataset consisting of just the observations with non-zero weights (i.e. the trial participants)
and assigning each observation a random number to assist in scrambling them for the plots, as well as
reducing the number of covariates included in the final data output*/

DATA _scrambleweights;
	SET _withweights;
	WHERE weight ^= 0 AND weight ^= .;
	CALL STREAMINIT (124);
	forscramble=RAND('UNIFORM');
	KEEP forscramble weight;
RUN;

PROC SORT DATA=_scrambleweights; /*Reordering the data by our random variable*/
	BY forscramble;
RUN;

DATA _forskyscraper;
	SET _scrambleweights;
	RETAIN place;
	IF _N_ = 1 THEN place = 0;
	ELSE place = place + weight;
	OUTPUT;
	KEEP weight place;
RUN;

ODS GRAPHICS / reset=index imagename="&FIGURENAME"; /*Assigning the image name*/

/*Now, we plot*/

PROC SGPLOT DATA=_forskyscraper;
	SERIES X=place Y=weight;
	YAXIS LABEL="Stabilized weights" MAX=20 MIN=0 LABELATTRS=(size=14) VALUEATTRS=(size=12);
	REFLINE 10 / axis=Y transparency = 0.5;
	XAXIS DISPLAY=(NOLABEL) LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;

%mend;

/*Executing the macro*/;

%Skyscraper(TRIALDATA=trialset,
		TARGETDATA=targset,
		CATPARAMS=female colon livermets WildKRAS,
		MODELPARAMS=age female colon livermets WildKRAS,
		FIGURENAME=Skyscraper);
