#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"


void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char *ret, tag[MAXTAGS][50];
  int errorFlag = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaBaryon");
  addr[nt] = &OmegaBaryon;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaDM_2ndSpecies");
  addr[nt] = &OmegaDM_2ndSpecies;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "HubbleParam");
  addr[nt] = &HubbleParam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ShapeGamma");
  addr[nt] = &ShapeGamma;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "PrimordialIndex");
  addr[nt] = &PrimordialIndex;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Box");
  addr[nt] = &Box;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "klimitB");
  addr[nt] = &klimitB;
  id[nt++] = FLOAT;


  strcpy(tag[nt], "B1mpc");
  addr[nt] = &B1mpc;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "boost");
  addr[nt] = &boost;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "nB");
  addr[nt] = &nB;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "coeff");
  addr[nt] = &coeff;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "kappa_b");
  addr[nt] = &kappa_b;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "kappa_dm");
  addr[nt] = &kappa_dm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "p_shape");
  addr[nt] = &p_shape;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "lambda_pmf_b");
  addr[nt] = &lambda_pmf_b;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "lambda_pmf_dm");
  addr[nt] = &lambda_pmf_dm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "lambda_lcdm_dm");
  addr[nt] = &lambda_lcdm_dm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "lambda_lcdm_b");
  addr[nt] = &lambda_lcdm_b;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "same_phase");
  addr[nt] = &same_phase;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Redshift");
  addr[nt] = &Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nmesh");
  addr[nt] = &Nmesh;
  id[nt++] = INT;

  strcpy(tag[nt], "Nsample");
  addr[nt] = &Nsample;
  id[nt++] = INT;

  strcpy(tag[nt], "GlassFile");
  addr[nt] = GlassFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputSpectrum");
  addr[nt] = FileWithInputSpectrum;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputSpectrumDM");
  addr[nt] = FileWithInputSpectrumDM;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputSpectrumB");
  addr[nt] = FileWithInputSpectrumB;
  id[nt++] = STRING;


  strcpy(tag[nt], "FileWithInputSpectrumDMpmf");
  addr[nt] = FileWithInputSpectrumDMpmf;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputSpectrumBpmf");
  addr[nt] = FileWithInputSpectrumBpmf;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithTransfer");
  addr[nt] = FileWithTransfer;
  id[nt++] = STRING;


  strcpy(tag[nt], "GlassTileFac");
  addr[nt] = &GlassTileFac;
  id[nt++] = INT;

  strcpy(tag[nt], "Seed");
  addr[nt] = &Seed;
  id[nt++] = INT;

  strcpy(tag[nt], "Seedb");
  addr[nt] = &Seedb;
  id[nt++] = INT;

  



  strcpy(tag[nt], "SphereMode");
  addr[nt] = &SphereMode;
  id[nt++] = INT;

  strcpy(tag[nt], "NumFilesWrittenInParallel");
  addr[nt] = &NumFilesWrittenInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileBase");
  addr[nt] = FileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "WhichSpectrum");
  addr[nt] = &WhichSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
  addr[nt] = &InputSpectrum_UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ReNormalizeInputSpectrum");
  addr[nt] = &ReNormalizeInputSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_On");
  addr[nt] = &WDM_On;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_Vtherm_On");
  addr[nt] = &WDM_Vtherm_On;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_PartMass_in_kev");
  addr[nt] = &WDM_PartMass_in_kev;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "NU_On");
  addr[nt] = &NU_On;
  id[nt++] = INT;

  strcpy(tag[nt], "NU_Vtherm_On");
  addr[nt] = &NU_Vtherm_On;
  id[nt++] = INT;

  strcpy(tag[nt], "NU_PartMass_in_ev");
  addr[nt] = &NU_PartMass_in_ev;
  id[nt++] = FLOAT;

  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  buf[0] = 0;
	  ret = fgets(buf, 200, fd);

	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy(addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
	      if(ThisTask == 0)
		fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname,
			buf1);
	      errorFlag = 1;
	    }
	}
      fclose(fd);

    }
  else
    {
      if(ThisTask == 0)
	fprintf(stdout, "Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	  if(ThisTask == 0)
	    fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

    const char* delimiter = "/";
    const char* extension = ".param";
    // Allocate a new buffer to store the extracted substring
    char extracted_substring[100]; // Adjust the size as needed
    // Find the last occurrence of the delimiter "/"
    char* last_occurrence = strrchr(fname, delimiter[0]);

    if (last_occurrence != NULL) {
        // Move the pointer one position ahead to skip the delimiter
        last_occurrence++;

        // Find the position of ".param" in the substring
        char* extension_position = strstr(last_occurrence, extension);

        if (extension_position != NULL) {
            // Calculate the length of the substring
            size_t length = extension_position - last_occurrence;

            // Copy the substring to the new buffer
            strncpy(extracted_substring, last_occurrence, length);
            extracted_substring[length] = '\0'; // Add null-terminator

            // printf("Original string: %s\n", fname);
            // printf("Substring after last '/' and before '.param': %s\n", extracted_substring);
        } else {
            printf("'.param' not found in the string after the last '/'.\n");
        }
    } else {
        printf("Delimiter '/' not found in the string.\n");
    }

  printf("\n output directory is ");
  strcat(OutputDir,extracted_substring);
  printf(OutputDir);
  printf("\n");
#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
