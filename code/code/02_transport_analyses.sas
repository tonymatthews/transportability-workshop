/*Running the transport analysis*/

/*Naming the library where the trial data are being kept*/


LIBNAME cohort "C:\Users\webst\OneDrive\Documents\SER seminar data\Trial";

/*Assigning the library where the target population data is held*/

LIBNAME targ "C:\Users\webst\OneDrive\Documents\SER seminar data\Target";

/*libname lwork slibref=work server=server;
libname lcohort slibref=cohort server=server;
libname ltarg slibref=targ server=server;*/


/*First, let's bring the trial data into the "work" directory*/

DATA trialset;
	SET cohort.analytictrialcohort;
RUN;

/*Repeating something similar for the target data*/

DATA targset;
	SET targ.analytictargcohort;
RUN;

/*Now, we want to obtain weights so that the trial resembles the target using the below macro*/;

%macro getweights(	TRIALDATA=,/*The name of the trial data set*/
					TARGETDATA=,/*The name of the target data set*/
					CATPARAMS=,/*The categorical variables for the model*/
					MODELPARAMS=,/*All the variables to be used in the model*/
					WEIGHTBYARM="Yes",/*Yes (reweight each arm to the target) or No (reweight the whole trial), defaults to Yes*/
					OUTPUTNAME=);/*The name of the output data set that includes weights*/

					
					

/*IF WEIGHTBYARM is "Yes," we need to split the trial data into two groups, one for each arm,
then combine each of them with the target data*/

%IF &WEIGHTBYARM="Yes" %THEN %DO;

DATA _arm1;
	SET &TRIALDATA;
	WHERE treatment=1;
RUN;

DATA _arm0;
	SET &TRIALDATA;
	WHERE treatment=0;
RUN;

DATA _concatarm1;
	SET _arm1 &TARGETDATA;
RUN;

DATA _concatarm0;
	SET _arm0 &TARGETDATA;
RUN;

/*Once we have our two concatenated data sets, we now run two logistic regressions including the
variables we want to standardize*/

PROC LOGISTIC DATA=_concatarm1 DESCENDING NOPRINT;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS ;
	OUTPUT OUT=_arm1preds pred=denom;
RUN;

PROC LOGISTIC DATA=_concatarm0 DESCENDING NOPRINT;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS ;
	OUTPUT OUT=_arm0preds pred=denom;
RUN;

/*To obtain the numerator of stabilized odds weights, we have to run two additional regressions
with an empty set of variables.*/

PROC LOGISTIC DATA=_arm1preds DESCENDING NOPRINT;
	MODEL intrial = ;
	OUTPUT OUT=_arm1preds pred=num;
RUN;

PROC LOGISTIC DATA=_arm0preds DESCENDING NOPRINT;
	MODEL intrial = ;
	OUTPUT OUT=_arm0preds pred=num;
RUN;

/*Next, we combine the data sets; note that those with intrial = 0 are present twice in this set*/

DATA _preweighted;
	SET _arm1preds _arm0preds;
RUN;

%END;

%ELSE %DO;

/*If WEIGHTBYARM is not "Yes", we instead simply concatenate the trial and target populations,
then run similar logistic regressions*/

DATA _concat;
	SET &TRIALDATA &TARGETDATA;
RUN;

PROC LOGISTIC DATA=_concat DESCENDING NOPRINT;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS ;
	OUTPUT OUT=_concatpred pred=denom;
RUN;

PROC LOGISTIC DATA=_concatpred DESCENDING NOPRINT;
	MODEL intrial = ;
	OUTPUT OUT=_preweighted pred=num;
RUN;

%END;

/*Now, we simply calculate stabilized or unstabilized inverse odds weights*/

DATA &OUTPUTNAME;
	SET _preweighted;
	denomodds = denom/(1-denom); /*First we convert the denominator from a probability to an odds*/
	numodds = num/(1-num); /*We do the same for the numerator of the stabilized weights*/
	IF intrial = 1 THEN DO; /*Assigning weights to trial participants*/
		iosw = 1/denomodds;
		siosw = numodds/denomodds;
	END;
	ELSE DO; /*Assigning zero weights to target participants*/
		iosw = 0;
		siosw = 0;
	END;
RUN;

%MEND;

%getweights(TRIALDATA=trialset, TARGETDATA=targset, CATPARAMS=female colon livermets WildKRAS, MODELPARAMS=age female colon livermets WildKRAS, WEIGHTBYARM="No",OUTPUTNAME=Outset);

/*Now that we have weights, we can use them in a weighted analysis to estimate a hazard ratio for the effect of treatment on progression-free survival*/

PROC PHREG DATA=Outset COVS; /*Specifying a sandwich estimator for the variance-note that to incorporate uncertainty in the target population, you need a bootstrap*/
	WHERE intrial = 1; /*Restricting to those in the trial; this is not necessary since target individuals have weights of 0 but is a useful reminder*/
	CLASS treatment (ref="0") / param=ref;
	MODEL pfsdycr*pfscr(0)=treatment;
	WEIGHT SIOSW;
RUN;

/*We can compare those results to the unweighted estimate of the hazard ratio*/

PROC PHREG DATA=Outset;
	WHERE intrial = 1; /*Restricting to those in the trial; this is not necessary since target individuals have weights of 0 but is a useful reminder*/
	CLASS treatment (ref="0") / param=ref;
	MODEL pfsdycr*pfscr(0)=treatment;
RUN;

/*If we wanted to run a bootstrap, we would do something like this*/

%macro bootstrapiosw(TRIALDATA=,
					TARGETDATA=,
					CATPARAMS=,
					MODELPARAMS=,
					WEIGHTBYARM="Yes",
					ITERATIONS=);

/*First, replicate the trial and target data for some number of iterations to get a sense of uncertainty in estimates */

PROC SURVEYSELECT DATA=&TRIALDATA OUT=_boottrial method=urs reps=&ITERATIONS seed=12345 samprate=1 OUTHITS;
RUN;

PROC SURVEYSELECT DATA=&TARGETDATA OUT=_boottarg method=urs reps=&ITERATIONS seed=135 samprate=1 OUTHITS;
RUN;

/*Next, add back in the original trial and target data as replicate 0*/

DATA _boottrial;
	SET &TRIALDATA _boottrial;
	IF replicate = . THEN replicate = 0;
	intrial = 1; /*Giving the study/trial participants a variable indicating their origin*/
RUN;

DATA _boottarg;
	SET &TARGETDATA _boottarg;
	IF replicate = . THEN replicate = 0;
	intrial = 0; /*Giving the target participants a variable indicating their origin*/
RUN;

/*Now, things split again depending on whether we selected by-arm weights or not*/


%IF &WEIGHTBYARM="Yes" %THEN %DO;

DATA _arm1;
	SET _boottrial;
	WHERE treatment=1;
RUN;

DATA _arm0;
	SET _boottrial;
	WHERE treatment=0;
RUN;

DATA _concatarm1;
	SET _arm1 _boottarg;
	BY replicate;
RUN;

DATA _concatarm0;
	SET _arm0 _boottarg;
	BY replicate;
RUN;

/*Once we have our two concatenated data sets, we now run two logistic regressions including the
variables we want to standardize*/

PROC LOGISTIC DATA=_concatarm1 DESCENDING NOPRINT;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS ;
	BY replicate;
	OUTPUT OUT=_arm1preds pred=denom;
RUN;

PROC LOGISTIC DATA=_concatarm0 DESCENDING NOPRINT;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS ;
	BY replicate;
	OUTPUT OUT=_arm0preds pred=denom;
RUN;

/*To obtain the numerator of stabilized odds weights, we have to run two additional regressions
with an empty set of variables.*/

PROC LOGISTIC DATA=_arm1preds DESCENDING NOPRINT;
	MODEL intrial = ;
	BY replicate;
	OUTPUT OUT=_arm1preds pred=num;
RUN;

PROC LOGISTIC DATA=_arm0preds DESCENDING NOPRINT;
	MODEL intrial = ;
	BY replicate;
	OUTPUT OUT=_arm0preds pred=num;
RUN;

/*Next, we combine the data sets; note that those with intrial = 0 are present twice in this set*/

DATA _preweighted;
	SET _arm1preds _arm0preds;
	BY replicate;
RUN;

%END;

%ELSE %DO;

/*If WEIGHTBYARM is not "Yes", we instead simply concatenate the trial and target populations,
then run similar logistic regressions*/

DATA _concat;
	SET _boottrial _boottarg;
	BY replicate;
RUN;

PROC LOGISTIC DATA=_concat DESCENDING NOPRINT;
	CLASS &CATPARAMS / param=ref;
	MODEL intrial = &MODELPARAMS ;
	BY replicate;
	OUTPUT OUT=_concatpred pred=denom;
RUN;

PROC LOGISTIC DATA=_concatpred DESCENDING NOPRINT;
	MODEL intrial = ;
	BY replicate;
	OUTPUT OUT=_preweighted pred=num;
RUN;

%END;

/*Now, we simply calculate stabilized or unstabilized inverse odds weights*/

DATA _withweights;
	SET _preweighted;
	denomodds = denom/(1-denom); /*First we convert the denominator from a probability to an odds*/
	numodds = num/(1-num); /*We do the same for the numerator of the stabilized weights*/
	IF intrial = 1 THEN DO; /*Assigning weights to trial participants*/
		iosw = 1/denomodds;
		siosw = numodds/denomodds;
	END;
	ELSE DO; /*Assigning zero weights to target participants*/
		iosw = 0;
		siosw = 0;
	END;
RUN;

/*And run our analyses*/

PROC PHREG DATA=_withweights COVS; /*Specifying a sandwich estimator for the variance-note that to incorporate uncertainty in the target population, you need a bootstrap*/
	WHERE intrial = 1; /*Restricting to those in the trial; this is not necessary since target individuals have weights of 0 but is a useful reminder*/
	BY replicate;
	CLASS treatment (ref="0") / param=ref;
	MODEL pfsdycr*pfscr(0)=treatment;
	WEIGHT SIOSW;
	ODS OUTPUT PARAMETERESTIMATES=parms;
RUN;

PROC MEANS DATA=parms;
	VAR Estimate;
	WHERE replicate ^= 0;
RUN;

%MEND;

%bootstrapiosw(TRIALDATA=trialset, TARGETDATA=targset, CATPARAMS=female colon livermets WildKRAS, MODELPARAMS=age female colon livermets WildKRAS, WEIGHTBYARM="No",ITERATIONS=250);
