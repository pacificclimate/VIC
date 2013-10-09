#include "WriteOutputBinary.h"

FILE* WriteOutputBinary::openFile(const char* string) {
  return open_file(string, "wb");
}

void WriteOutputBinary::write_data(out_data_file_struct* out_data_files,
    out_data_struct* out_data, const dmy_struct* dmy, int dt,
    const ProgramState* state) {
  // Initialize pointers
  char               *tmp_cptr = (char *)calloc(N_OUTVAR_TYPES*state->options.Nlayer*state->options.SNOW_BAND,sizeof(char));
  short int          *tmp_siptr = (short int *)calloc(N_OUTVAR_TYPES*state->options.Nlayer*state->options.SNOW_BAND,sizeof(short int));
  unsigned short int *tmp_usiptr = (unsigned short int *)calloc(N_OUTVAR_TYPES*state->options.Nlayer*state->options.SNOW_BAND,sizeof(unsigned short int));
  int                *tmp_iptr = (int *)calloc(N_OUTVAR_TYPES*state->options.Nlayer*state->options.SNOW_BAND,sizeof(int));
  float              *tmp_fptr = (float *)calloc(N_OUTVAR_TYPES*state->options.Nlayer*state->options.SNOW_BAND,sizeof(float));
  double             *tmp_dptr = (double *)calloc(N_OUTVAR_TYPES*state->options.Nlayer*state->options.SNOW_BAND,sizeof(double));

  // Time
  tmp_iptr[0] = dmy->year;
  tmp_iptr[1] = dmy->month;
  tmp_iptr[2] = dmy->day;
  tmp_iptr[3] = dmy->hour;

  // Loop over output files
  for (int file_idx = 0; file_idx < state->options.Noutfiles; file_idx++) {

#if !OUTPUT_FORCE
    // Write the date
    if (dt < 24) {
      // Write year, month, day, and hour
      fwrite(tmp_iptr, sizeof(int), 4, out_data_files[file_idx].fh);
    }
    else {
      // Only write year, month, and day
      fwrite(tmp_iptr, sizeof(int), 3, out_data_files[file_idx].fh);
    }
#endif

    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < out_data_files[file_idx].nvars; var_idx++) {
      // Loop over this variable's elements
      int ptr_idx = 0;
      if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_CHAR) {
        for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          tmp_cptr[ptr_idx++] = (char)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
        }
        fwrite(tmp_cptr, sizeof(char), ptr_idx, out_data_files[file_idx].fh);
      }
      else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_SINT) {
        for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          tmp_siptr[ptr_idx++] = (short int)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
        }
        fwrite(tmp_siptr, sizeof(short int), ptr_idx, out_data_files[file_idx].fh);
      }
      else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_USINT) {
        for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          tmp_usiptr[ptr_idx++] = (unsigned short int)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
        }
        fwrite(tmp_usiptr, sizeof(unsigned short int), ptr_idx, out_data_files[file_idx].fh);
      }
      else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_INT) {
        for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          tmp_iptr[ptr_idx++] = (int)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
        }
        fwrite(tmp_iptr, sizeof(int), ptr_idx, out_data_files[file_idx].fh);
      }
      else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_FLOAT) {
        for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          tmp_fptr[ptr_idx++] = (float)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
        }
        fwrite(tmp_fptr, sizeof(float), ptr_idx, out_data_files[file_idx].fh);
      }
      else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_DOUBLE) {
        for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          tmp_dptr[ptr_idx++] = (double)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
        }
        fwrite(tmp_dptr, sizeof(double), ptr_idx, out_data_files[file_idx].fh);
      }
    }

  }

  // Free the arrays
  free((char *)tmp_cptr);
  free((char *)tmp_siptr);
  free((char *)tmp_usiptr);
  free((char *)tmp_iptr);
  free((char *)tmp_fptr);
  free((char *)tmp_dptr);
}

/* Binary header format:

 Data        Stored As           Comment

 Identifier  (unsigned short)*4  0xFFFF, repeated 4 times
 Nbytes      (unsigned short)*1  Number of bytes in the header,
                                 INCLUDING THE IDENTIFIER

 Part 1: Global Attributes
 Nbytes1     (unsigned short)*1  Number of bytes in part 1
 nrecs       (int)*1             Number of records in the file
 dt          (int)*1             Output time step length in hours
 startyear   (int)*1             Year of first record
 startmonth  (int)*1             Month of first record
 startday    (int)*1             Day of first record
 starthour   (int)*1             Hour of first record
 ALMA_OUTPUT (char)*1            0 = standard VIC units; 1 = ALMA units
 Nvars       (char)*1            Number of variables in the file, including date fields

 Part 2: Variables
 Nbytes2     (unsigned short)*1  Number of bytes in part 2
 For each variable, the following fields: { len varname type mult }
   len       (char)*1            Number of characters in varname
   varname   (char)*len          Variable name
   type      (char)*1            Code identifying variable type
   mult      (float)*1           Multiplier for variable*/
void WriteOutputBinary::write_header(out_data_file_struct* out_data_files,
    out_data_struct* out_data, const dmy_struct* dmy,
    const ProgramState* state) {
  unsigned short      Identifier;
  unsigned short      Nbytes;
  unsigned short      Nbytes1;
  unsigned short      Nbytes2;
  char                tmp_ALMA_OUTPUT;
  char                Nvars;
  char                tmp_len;
  char                tmp_type;
  float               tmp_mult;

  prepareDataForWriting(out_data);

  if (state->options.ALMA_OUTPUT)
    tmp_ALMA_OUTPUT = 1;
  else
    tmp_ALMA_OUTPUT = 0;

  char tmp_str[256];

  // Identifier
  Identifier = 0xFFFF;

  // Loop over output files
  for (int file_idx = 0; file_idx < state->options.Noutfiles; file_idx++) {

    // ***** Compute the number of bytes in part 1 *****

    // 1 instance of Nbytes1
    Nbytes1 = sizeof(unsigned short);

    // nrecs
    Nbytes1 += sizeof(int);

    // dt
    Nbytes1 += sizeof(int);

    // start date (year, month, day, hour)
    Nbytes1 += 4*sizeof(int);

    // ALMA_OUTPUT
    Nbytes1 += sizeof(char);

    // Nvars
    Nbytes1 += sizeof(char);

    // ***** Compute the number of bytes in part 2 *****

    // 1 instance of Nbytes2
    Nbytes2 = sizeof(unsigned short);

#if !OUTPUT_FORCE
    // Date fields
    Nbytes2 += sizeof(char) + 4*sizeof(char) + sizeof(char) + sizeof(float); // year
    Nbytes2 += sizeof(char) + 5*sizeof(char) + sizeof(char) + sizeof(float); // month
    Nbytes2 += sizeof(char) + 3*sizeof(char) + sizeof(char) + sizeof(float); // day
    if (state->global_param.out_dt < 24)
      Nbytes2 += sizeof(char) + 4*sizeof(char) + sizeof(char) + sizeof(float); // hour
#endif

    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < out_data_files[file_idx].nvars; var_idx++) {
      // Loop over this variable's elements
      for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
        if (out_data[out_data_files[file_idx].varid[var_idx]].nelem > 1)
          sprintf(tmp_str, "%s_%d", out_data[out_data_files[file_idx].varid[var_idx]].varname, elem_idx);
        else
          strcpy(tmp_str, out_data[out_data_files[file_idx].varid[var_idx]].varname);
        Nbytes2 += sizeof(char) + strlen(tmp_str)*sizeof(char) + sizeof(char) + sizeof(float);
      }
    }

    // ***** Compute the total number of bytes in the header *****

    // 4 instances of Identifier, plus 1 instance of Nbytes, plus number of bytes in parts 1 and 2
    Nbytes = 4*sizeof(unsigned short) + sizeof(unsigned short) + Nbytes1 + Nbytes2;

    // ***** Write the header *****

    // 4 instances of Identifier
    for (int i=0; i<4; i++)
      fwrite(&Identifier, sizeof(unsigned short), 1, out_data_files[file_idx].fh);

    // Nbytes
    fwrite(&Nbytes, sizeof(unsigned short), 1, out_data_files[file_idx].fh);

    // Nbytes1
    fwrite(&Nbytes1, sizeof(unsigned short), 1, out_data_files[file_idx].fh);

    // nrecs
    fwrite(&(state->global_param.nrecs), sizeof(int), 1, out_data_files[file_idx].fh);

    // dt
    fwrite(&(state->global_param.out_dt), sizeof(int), 1, out_data_files[file_idx].fh);

    // start date (year, month, day, hour)
    fwrite(&(dmy->year), sizeof(int), 1, out_data_files[file_idx].fh);
    fwrite(&(dmy->month), sizeof(int), 1, out_data_files[file_idx].fh);
    fwrite(&(dmy->day), sizeof(int), 1, out_data_files[file_idx].fh);
    fwrite(&(dmy->hour), sizeof(int), 1, out_data_files[file_idx].fh);

    // ALMA_OUTPUT
    fwrite(&tmp_ALMA_OUTPUT, sizeof(char), 1, out_data_files[file_idx].fh);

    // Nvars
    Nvars = out_data_files[file_idx].nvars;
#if !OUTPUT_FORCE
    if (state->global_param.out_dt < 24)
      Nvars += 4;
    else
      Nvars += 3;
#endif
    fwrite(&Nvars, sizeof(char), 1, out_data_files[file_idx].fh);

    // Nbytes2
    fwrite(&Nbytes2, sizeof(unsigned short), 1, out_data_files[file_idx].fh);

#if !OUTPUT_FORCE
    // Date fields
    tmp_type = OUT_TYPE_INT;
    tmp_mult = 1.;

    // year
    strcpy(tmp_str,"YEAR");
    tmp_len = strlen(tmp_str);
    fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
    fwrite(tmp_str, sizeof(char), tmp_len, out_data_files[file_idx].fh);
    fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
    fwrite(&tmp_mult, sizeof(float), 1, out_data_files[file_idx].fh);

    // month
    strcpy(tmp_str,"MONTH");
    tmp_len = strlen(tmp_str);
    fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
    fwrite(tmp_str, sizeof(char), tmp_len, out_data_files[file_idx].fh);
    fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
    fwrite(&tmp_mult, sizeof(float), 1, out_data_files[file_idx].fh);

    // day
    strcpy(tmp_str,"DAY");
    tmp_len = strlen(tmp_str);
    fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
    fwrite(tmp_str, sizeof(char), tmp_len, out_data_files[file_idx].fh);
    fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
    fwrite(&tmp_mult, sizeof(float), 1, out_data_files[file_idx].fh);

    if (state->global_param.out_dt < 24) {
      // hour
      strcpy(tmp_str,"HOUR");
      tmp_len = strlen(tmp_str);
      fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
      fwrite(tmp_str, sizeof(char), tmp_len, out_data_files[file_idx].fh);
      fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
      fwrite(&tmp_mult, sizeof(float), 1, out_data_files[file_idx].fh);
    }

#endif

    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < out_data_files[file_idx].nvars; var_idx++) {
      // Loop over this variable's elements
      for (int elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
        if (out_data[out_data_files[file_idx].varid[var_idx]].nelem > 1)
          sprintf(tmp_str, "%s_%d", out_data[out_data_files[file_idx].varid[var_idx]].varname, elem_idx);
        else
          strcpy(tmp_str, out_data[out_data_files[file_idx].varid[var_idx]].varname);
        tmp_len = strlen(tmp_str);
        fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
        fwrite(tmp_str, sizeof(char), tmp_len, out_data_files[file_idx].fh);
        tmp_type = out_data[out_data_files[file_idx].varid[var_idx]].type;
        fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
        tmp_mult = out_data[out_data_files[file_idx].varid[var_idx]].mult;
        fwrite(&tmp_mult, sizeof(float), 1, out_data_files[file_idx].fh);
      }
    }
  }
}

const char* WriteOutputBinary::getDescriptionOfOutputType() {
  return "BINARY";
}

void WriteOutputBinary::prepareDataForWriting(out_data_struct* out_data) {
  for (int v=0; v<N_OUTVAR_TYPES; v++) {
    for (int i=0; i<out_data[v].nelem; i++) {
      out_data[v].aggdata[i] *= out_data[v].mult;
    }
  }
}



