/*Visualization #2: VITT plot*/

/*Naming the library where the trial data are being kept*/

LIBNAME cohort "C:\Users\webst\OneDrive\Documents\SER seminar data\Trial";

/*Naming the library where the target data are being kept*/

LIBNAME targ "C:\Users\webst\OneDrive\Documents\SER seminar data\Target";

ods listing gpath= "C:\Users\webst\OneDrive\Documents\SER seminar data\Visualizations"; /*Assigning folder where figures will be saved*/

/*Creating some data sets to designate colors and formatting*/

DATA customattrs;
length markercolor $ 9 value $ 20;
input ID $ value $ markercolor $;
datalines;
myid  AGE red
myid  WildKRAS yellow
myid  colon green
myid  female blue
myid  livermets purple
;
RUN;

DATA customattrs2;
	SET customattrs;
	IF value="AGE" THEN value = "Age/20";
	IF value="WildKRAS" THEN value="Wild-type KRAS";
	IF value ="colon" THEN value="Colon cancer";
	IF value="female" THEN value="Female sex";
	IF value="livermets" THEN value="Liver metastases";
RUN;

DATA _delta;
	LENGTH parameterlegend $ 20;
	xORsamp = 1000;
	xORout = 1000;
	parameterlegend = "Age/20";
	OUTPUT;
	parameterlegend = "Liver metastases";
	OUTPUT;
	parameterlegend = "Wild-type KRAS";
	OUTPUT;
	parameterlegend = "Female sex";
	OUTPUT;
	parameterlegend = "Colon cancer";
	OUTPUT;
RUN;

/*First, let's bring the trial data into the "work" directory, and recode to a binary 1/0 outcome for one-year progression free survival; let's
also rework age to have units of about one standard deviation of the data, or 20 years*/

DATA trialset;
	SET cohort.analytictrialcohort;
	IF pfsdycr >= 365 THEN outcome = 0; /*Those with >365 days of survival didn't have the outcome*/
	ELSE IF pfscr = 1 THEN outcome = 1; /*Those with <365 days of survival who had the outcome had the outcome*/
	ELSE outcome = .; /*Other individuals were censored for the outcome*/
	age=age/20;
RUN;

/*Repeating something similar for the target data*/

DATA targset;
	SET targ.analytictargcohort;
	age=age/20;
RUN;


ods listing gpath= "/local/users/mawc/PCORI/SER Seminar/Visualizations"; /*Assigning folder where figures will be saved*/


/*Next, getting the quantities we need for the chart*/;

%macro VITTplot(	TRIALDATA=, /*Name of trial/study population file*/
				TARGETDATA=, /*Name of target population file*/
				MODELPARAMS=,/*All the variables to be used in the model*/
				TREATMENTINOUTCOME="Yes", /*Specify whether treatment will be included in the outcome model or not-may want to examine both*/
				FIGURENAME=);

/*First, replicate the trial and target data 500 times to get a sense of uncertainty in estimates */

PROC SURVEYSELECT DATA=&TRIALDATA OUT=_boottrial method=urs reps=500 seed=12345 samprate=1;
RUN;

PROC SURVEYSELECT DATA=&TARGETDATA OUT=_boottarg method=urs reps=500 seed=1235 samprate=1;
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

/*Then, use the trial data to estimate multivariable odds ratios between variables and the outcome*/


ods graphics off;

ods exclude all;

ods noresults;

%IF &TREATMENTINOUTCOME="Yes" %THEN %DO;

PROC GENMOD DATA=_boottrial DESCENDING;
	BY replicate;
	MODEL outcome = &MODELPARAMS  treatment / link=logit dist=binomial;
	ODS OUTPUT PARAMETERESTIMATES=_outparams; /*Saving the results as "_outparams"*/
RUN;

%END;

%ELSE %DO;

PROC GENMOD DATA=_boottrial DESCENDING;
	BY replicate;
	MODEL outcome = &MODELPARAMS  / link=logit dist=binomial;
	ODS OUTPUT PARAMETERESTIMATES=_outparams; /*Saving the results as "_outparams"*/
RUN;

%END;

/*Next, we concatenate the study/trial and target data for estimating the sampling OR*/

DATA _concat;
	SET _boottrial _boottarg;
	BY replicate;
RUN;

/*And we run another logistic regression, this time on the intrial variable*/

PROC GENMOD DATA=_concat DESCENDING;
	BY replicate;
	MODEL intrial = &MODELPARAMS / link=logit dist=binomial;
	ODS OUTPUT PARAMETERESTIMATES=_sampparams;
RUN;


ods graphics on;

ods exclude none;

ods results;

PROC SQL;
	create table _prechart as
		select a.parameter, a.replicate, a.estimate as logORout, b.estimate as logORsamp from
			_outparams as a inner join _sampparams as b
				on a.parameter = b.parameter and a.replicate = b.replicate and a.parameter ^= "Intercept" and a.parameter ^= "Scale";
QUIT;

DATA _forchart;
	SET _prechart;
	LENGTH parameter2 $ 25.;

	ORout = exp(logORout); /*Exponentiating the estimated variable-outcome ORs*/
	ORsamp = exp(logORsamp); /*Exponentiating the estimated variable-sampling ORs*/

	IF replicate = 0 THEN DO; /*Separating the "base" odds ratios from the bootstrapped ones*/
		coreORout = ORout;
		coreORsamp = ORsamp;
	END;

	ELSE DO; /*Continuing this process*/
		repORout = ORout;
		repORsamp = ORsamp;
	END;	
	
	IF parameter="AGE" THEN parameter2 = "Age/20";
	IF parameter="WildKRAS" THEN parameter2="Wild-type KRAS";
	IF parameter ="colon" THEN parameter2="Colon cancer";
	IF parameter="female" THEN parameter2="Female sex";
	IF parameter="livermets" THEN parameter2="Liver metastases";

	parameterlegend="Age/20";
RUN;

DATA _finalchart;
	SET _delta _forchart;
RUN;

ODS GRAPHICS / reset=index imagename="&FIGURENAME"; /*Assigning the image name*/

PROC SGPLOT DATA=_finalchart dattrmap=customattrs2 aspect=1;
	SCATTER X=repORsamp Y=repORout / group=parameter2 markerattrs=(symbol = circlefilled) transparency = 0.8 attrid=myid;
	SCATTER X=coreORsamp Y=coreORout / markerattrs=(symbol=diamondfilled color=black);
	SCATTER X=xORsamp Y=xORout / group=parameterlegend markerattrs=(symbol=circlefilled size=10) attrid=myid name="label";
	REFLINE 1 / AXIS=X;
	REFLINE 1 / AXIS=Y;
	YAXIS LABEL="Outcome odds ratios" MAX=3.0 MIN=0.33 TYPE=LOG LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABEL="'Sampling' odds ratios" Max=3.0 MIN=0.33 TYPE=LOG LABELATTRS=(size=14) VALUEATTRS=(size=12);
	KEYLEGEND "label" / TITLE="Parameters";
RUN;

%mend;

/*Executing the macro*/;

%VITTPLOT(TRIALDATA=trialset,
		TARGETDATA=targset,
		MODELPARAMS=age female colon livermets WildKRAS,
		TREATMENTINOUTCOME=Yes,
		FIGURENAME=VITT_plot);

PROC SGPLOT DATA=_finalchart dattrmap=customattrs2 aspect=1;
	SCATTER X=repORsamp Y=repORout / group=parameter2 markerattrs=(symbol = circlefilled) transparency = 0.8 attrid=myid;
	SCATTER X=coreORsamp Y=coreORout / markerattrs=(symbol=diamondfilled color=black);
	SCATTER X=xORsamp Y=xORout / group=parameterlegend markerattrs=(symbol=circlefilled size=10) attrid=myid name="label";
	REFLINE 1 / AXIS=X;
	REFLINE 1 / AXIS=Y;
	YAXIS LABEL="Outcome odds ratios" MAX=8 MIN=0.125 TYPE=LOG LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABEL="'Sampling' odds ratios" Max=8 MIN=0.125 TYPE=LOG LABELATTRS=(size=14) VALUEATTRS=(size=12);
	KEYLEGEND "label" / TITLE="Parameters";
RUN;

PROC SGPLOT DATA=_finalchart dattrmap=customattrs2 aspect=1;
	WHERE parameter = "WildKRAS" AND replicate ^= 0;
	SCATTER X=repORsamp Y=repORout / group=parameter2 markerattrs=(symbol = circlefilled) transparency = 0.8 attrid=myid;
	SCATTER X=coreORsamp Y=coreORout / markerattrs=(symbol=diamondfilled color=black);
	SCATTER X=xORsamp Y=xORout / group=parameterlegend markerattrs=(symbol=circlefilled size=10) attrid=myid name="label";
	REFLINE 1 / AXIS=X;
	REFLINE 1 / AXIS=Y;
	YAXIS LABEL="Outcome odds ratios" MAX=8 MIN=0.125 TYPE=LOG LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABEL="'Sampling' odds ratios" Max=8 MIN=0.125 TYPE=LOG LABELATTRS=(size=14) VALUEATTRS=(size=12);
	KEYLEGEND "label" / TITLE="Parameters";
RUN;
