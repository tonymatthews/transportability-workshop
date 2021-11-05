/*Trimming down and reorganizing the trial and target population data*/


OPTIONS NOFMTERR;

/*Assigning the library where the cohort data is held*/

LIBNAME cohort "C:\Users\webst\OneDrive\Documents\SER seminar data\Trial";

/*Assigning the library where the target population data is held*/

LIBNAME targ "C:\Users\webst\OneDrive\Documents\SER seminar data\Target";


/*First, we use a PROC SQL statement to join the covariate and KRAS trial data sets*/

PROC SQL;
	create table trialcohortv1 as
		select a.*, b.BMMTR1 as KRAS from
			cohort.adsl_pds2019 as a left join cohort.biomark_pds2019 as b
				on a.subjid = b.subjid;
QUIT;

/*Next, we turn that data into our analytic trial data set, trimming some extraneous covariates and renaming others*/

DATA cohort.analytictrialcohort;
	SET trialcohortv1;

	intrial = 1; /*This is our trial indicator variable*/


	/*Now to recode some key covariates to numerical quantities to make some figures easier to make*/

	IF KRAS = " " OR KRAS = "Failure" THEN WildKRAS = .; /*First, wild-type KRAS*/
	ELSE IF KRAS = "Wild-type" THEN WildKRAS = 1;
	ELSE WildKRAS = 0;


	IF trt="FOLFOX alone" THEN treatment = 0; /*Next, treatment*/
	ELSE treatment = 1;

	IF diagtype = "Colon" THEN colon = 1; /*Next, colon cancer site (colon vs rectum)*/
	ELSE colon = 0;

	IF sex = "Male" THEN female = 0; /*Next, sex*/
	ELSE female = 1;

	IF livermet = "Y" THEN livermets = 1; /*And liver metastases*/
	ELSE livermets = 0;

	/*Next, we output those with no missing values for our new covariates*/

	IF WildKRAS ^= . AND female ^= . AND colon ^= . AND livermets ^= . AND age ^= . THEN OUTPUT;

	/*Now we just keep our new variables, as well as the progression-free survival data*/

	KEEP intrial treatment WildKRAS female colon livermets  age pfsdycr pfscr;
RUN;

/*Next, renaming the target cohort we're working with*/

DATA targ.analytictargcohort;
	SET targ.basepop;

	/*Assinging a 0 value to the intrial variable*/

	intrial = 0;
	KEEP intrial female colon livermets age WildKRAS;
RUN;
