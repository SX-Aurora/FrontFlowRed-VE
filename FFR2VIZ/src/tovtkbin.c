/**
 * @author    Takashi Matsushita
 * @date      $Date: 2005/10/13 09:14:49 $
 * @version   $Revision: 1.4 $
 */
/*
 * Copyright: (c) 2005 Takashi Matsushita, All rights reserved.
 * License:   none
 * Created:   30 September 2005
 * Credits:   none
 */

/** @todo     nope */
/** @warnings nope */

/**
 * ffr to vtk converter (binary)
 */
/*=================================================================*
 * declarations
 *=================================================================*/
/*-----------------------------------------------------------------* 
 * constants
 *-----------------------------------------------------------------*/
#define VTKDATA_V2 "# vtk DataFile Version 2.0"
#define VTKDATA_ASCII "ASCII"
#define VTKDATA_BINARY "BINARY"
#define VTKDATA_UNSTRUCTURED_GRID "DATASET UNSTRUCTURED_GRID"
#define VTKDATA_POINT "POINTS"
#define VTKDATA_CELL "CELLS"
#define VTKDATA_CELL_TYPE "CELL_TYPES"
#define VTKDATA_CELL_DATA "CELL_DATA"
#define VTKDATA_VECTOR "VECTORS"
#define VTKDATA_SCALAR "SCALARS"
#define VTKDATA_LOOKUP_TABLE "LOOKUP_TABLE"
#define VTKDATA_POINT_DATA "POINT_DATA"

#define VTKDATA_TETRA 10
#define VTKDATA_HEXAHEDRON 12
#define VTKDATA_WEDGE 13
#define VTKDATA_PYRAMID 14

#define NODES_PER_TETRA 4
#define NODES_PER_HEXA 8
#define NODES_PER_WEDGE 6
#define NODES_PER_PYRAMID 5
#define MAX_NODES_PER_CELL 8

#define BINARY
#if 0 /* for debug */
#undef BINARY
#define VERBOSE
#endif

/*-----------------------------------------------------------------* 
 * headers
 *-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/*-----------------------------------------------------------------* 
 * file scope variables
 *-----------------------------------------------------------------*/
static const char rcsid[] = "$Id: tovtkbin.c,v 1.4 2005/10/13 09:14:49 takashi Exp takashi $";

/*-----------------------------------------------------------------* 
 * macros
 *-----------------------------------------------------------------*/
#define RSTRIP(_str) \
  do { \
    size_t _len = strlen(_str); \
    while (_len--) { \
      if (!isspace(_str[_len])) \
        break; \
      _str[_len] = '\0'; \
    } \
  } while (0)

#define SYS_EXIT(_name) \
  do { \
    fprintf(stderr, "%s:%s\n", __FILE__, __LINE__); \
    perror(_name); \
    exit(EXIT_FAILURE); \
  } while (0)

#define BSWAP2_(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define BSWAP4_(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | \
                    (((x) >> 8) & 0x0000ff00) | ((x) >> 24) )

#define BSWAP2(x) (*(unsigned short *)&(x) = BSWAP2_(*(unsigned short *)&(x)))
#define BSWAP4(x) (*(unsigned *)&(x) = BSWAP4_(*(unsigned *)&(x)))



/*=================================================================*
 * implementation
 *=================================================================*/
/*-----------------------------------------------------------------* 
 * file scope function
 *-----------------------------------------------------------------*/
/**
 * checks endian of this machine
 *
 * @note    none
 * @param   none
 * @return  one if little endian, 0 if big endian
 */
static int
is_little_endian(void)
{
  short int word = 0x0001;
  char *byte = (char *)&word;

  return (byte[0] ? 1 : 0);
}

/**
 * opens a file to write
 *
 * @note    none
 * @param   name [in] name of a file to open
 * @param   mode [in] file open mode
 * @return  pointer to opend file stream
 *          call exit(EXIT_FAILURE) if error occured
 */
static FILE *
vtk_open(const char *name, const char *mode)
{
  FILE *file;

  if ((file = fopen(name, mode)) == NULL) {
    SYS_EXIT("fopen");
  }

  return file;
}

/**
 * closes a file
 *
 * @note    none
 * @param   file [in] pointer to file stream
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
static void
vtk_close(FILE *file)
{
  if (fclose(file) == EOF) {
    SYS_EXIT("fclose");
  }

  return;
}

/**
 * malloc wrapper
 *
 * @note    none
 * @param   len [in] length in byte to allocate
 * @return  pointer to an allocated memory region
 *          call exit(EXIT_FAILURE) if error occured
 */
static void *
vtk_malloc(size_t len)
{
  void *buf;

  if ((buf = malloc(len)) == NULL) {
    SYS_EXIT("malloc");
  }
  memset(buf, 0, len);

  return buf;
}



/*-----------------------------------------------------------------* 
 * public function
 *-----------------------------------------------------------------*/
/**
 * writes vtk header part
 *
 * @note    none
 * @param   name [in] output file name
 * @param   len [in] length of name string
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
void
#ifdef WIN32
WRITE_VTK_HEADER_BIN(const char *name, int len)
#else
write_vtk_header_bin_(const char *name, int len)
#endif
{
#undef FNAME
#define FNAME "write_vtk_header_bin_"

  FILE *file;
  char *buf;
  
  buf = vtk_malloc(len+1);
  strncpy(buf, name, len);
  RSTRIP(buf);

#ifdef VERBOSE
  fprintf(stdout, "dbg> %s: opening '%s'\n", FNAME, buf);
#endif

  file = vtk_open(buf, "wb");
  if (fprintf(file,
              "%s\ncomment\n%s\n\n%s\n",
              VTKDATA_V2,
#ifdef BINARY
              VTKDATA_BINARY,
#else
              VTKDATA_ASCII,
#endif
              VTKDATA_UNSTRUCTURED_GRID
              ) < 0) {
    SYS_EXIT("fprintf");
  }

  vtk_close(file);
  free(buf);

  return;
}

/**
 * writes coordinate points of nodes
 *
 * @note    none
 * @param   name [in] output file name
 * @param   length [in] pointer to length of array
 * @param   array [in] pointer to coordinate points
 * @param   len [in] length of name string
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
void
#ifdef WIN32
WRITE_VTK_POINTS_BIN(const char *name, int *length, double *array, int len)
#else
write_vtk_points_bin_(const char *name, int *length, double *array, int len)
#endif
{
#undef FNAME
#define FNAME "write_vtk_points_bin_"

  FILE *file;
  char *buf;
  float real;
  int ii;
  int isLittle = is_little_endian();
  
  buf = vtk_malloc(len+1);
  strncpy(buf, name, len);
  RSTRIP(buf);

#ifdef VERBOSE
  fprintf(stdout, "dbg> %s: opening '%s'\n", FNAME, buf);
#endif

  file = vtk_open(buf, "a+b");

  if (fprintf(file, "%s %d float\n", VTKDATA_POINT, *length) < 0) {
    SYS_EXIT("fprintf");
  }

#ifdef BINARY
  for (ii = 0; ii < *length * 3; ii++) {
    real = (float)array[ii];
    if (isLittle) {
      BSWAP4(real);
    }
    if (fwrite(&real, sizeof(real), 1, file) != 1) {
      SYS_EXIT("fwrite");
    }
  }

#else /* for debug */
  for (ii = 0; ii < *length; ii++) {
    if (fprintf(file,
                "%E %E %E\n",
                array[0 + 3*ii],
                array[1 + 3*ii],
                array[2 + 3*ii]) < 0) {
      SYS_EXIT("fprintf");
    }
  }
#endif

  vtk_close(file);
  free(buf);

  return;
}

/**
 * writes cell data
 *
 * @note    none
 * @param   name [in] output file name
 * @param   length [in] pointer to length of array
 * @param   type [in] pointer to cell type array
 * @param   cell [in] pointer to node connectivity array of each cell
 * @param   len [in] length of name string
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
void
#ifdef WIN32
WRITE_VTK_CELLS_BIN(const char *name,
                     int *length,
                     int *type,
                     int *cell,
                     int len)
#else
write_vtk_cells_bin_(const char *name,
                     int *length,
                     int *type,
                     int *cell,
                     int len)
#endif
{
#undef FNAME
#define FNAME "write_vtk_cells_bin_"

  FILE *file;
  char *buf;
  int ii, jj, count;
  int num_tetra = 0, num_hexa = 0, num_wedge = 0, num_pyramid = 0;
  int num_tokens;
  int nodes_map[MAX_NODES_PER_CELL];
  int swap;
  int isLittle = is_little_endian();
  
  fprintf(stdout, "dbg> %s: %d\n", FNAME, *length);

  /* count number of cells */
  for (ii = 0; ii < *length; ii++) {
    switch (type[ii]) {
      case 1:
        num_tetra++;
        break;
      case 2:
        num_hexa++;
        break;
      case 3:
        num_wedge++;
        break;
      case 4:
        num_pyramid++;
        break;
      default:
        fprintf(stderr, "err> %s:unknown type:%d\n", FNAME, type[ii]);
        exit(EXIT_FAILURE);
    }
  }

  fprintf(stdout, "dbg> %s: # of tetra cells = %d\n", FNAME, num_tetra);
  fprintf(stdout, "dbg> %s: # of hexa cells = %d\n", FNAME, num_hexa);
  fprintf(stdout, "dbg> %s: # of wedge cells = %d\n", FNAME, num_wedge);
  fprintf(stdout, "dbg> %s: # of pyramid cells = %d\n", FNAME, num_pyramid);

  num_tokens = num_tetra * (NODES_PER_TETRA + 1) +
               num_hexa * (NODES_PER_HEXA + 1) +
               num_wedge * (NODES_PER_WEDGE + 1) +
               num_pyramid * (NODES_PER_PYRAMID + 1);

  buf = vtk_malloc(len+1);
  strncpy(buf, name, len);
  RSTRIP(buf);

#ifdef VERBOSE
  fprintf(stdout, "dbg> %s: opening '%s'\n", FNAME, buf);
#endif

  file = vtk_open(buf, "a+b");

  if (fprintf(file, "\n%s %d %d\n", VTKDATA_CELL, *length, num_tokens) < 0) {
    SYS_EXIT("fprintf");
  }

  /* dump node connectivity */
  for (ii = 0; ii < *length; ii++) {
    for (jj = 0; jj < MAX_NODES_PER_CELL; jj++) {
      nodes_map[jj] = jj;
    }

    switch (type[ii]) {
      case 1:
        count = NODES_PER_TETRA;
        break;

      case 2:
        count = NODES_PER_HEXA;
        break;

      case 3:
        nodes_map[0] = 0;
        nodes_map[1] = 2;
        nodes_map[2] = 1;
        nodes_map[3] = 3;
        nodes_map[4] = 5;
        nodes_map[5] = 4;

        count = NODES_PER_WEDGE;
        break;

      case 4:
        count = NODES_PER_PYRAMID;
        break;

      default:
        fprintf(stderr, "err> %s:unknown type:%d\n", FNAME, type[ii]);
        exit(EXIT_FAILURE);
    }

#ifdef BINARY
    swap = count;
    if (isLittle) {
      BSWAP4(swap);
    }

    if (fwrite(&swap, sizeof(swap), 1, file) != 1) {
      SYS_EXIT("fwrite");
    }

    for (jj = 0; jj < count; jj++) {
      /* vtk has zero node offset */
      swap = cell[ii*MAX_NODES_PER_CELL + nodes_map[jj]] - 1;
      if (isLittle) {
        BSWAP4(swap);
      }
      if (fwrite(&swap, sizeof(swap), 1, file) != 1) {
        SYS_EXIT("fwrite");
      }
    }

#else
    if (fprintf(file, "%d ", count) < 0) {
      SYS_EXIT("fprintf");
    }

    for (jj = 0; jj < count; jj++) {
      /* vtk has zero node offset */
      if (fprintf(file,
                  "%d ",
                  cell[ii*MAX_NODES_PER_CELL + nodes_map[jj]] - 1) < 0) {
        SYS_EXIT("fprintf");
      }
    }

    if (putc('\n', file) == EOF) {
      SYS_EXIT("putc");
    }
#endif
  }

  /* dump cell type */
  if (fprintf(file, "\n%s %d\n", VTKDATA_CELL_TYPE, *length) < 0) {
    SYS_EXIT("fprintf");
  }

  for (ii = 0; ii < *length; ii++) {
    switch (type[ii]) {
      case 1:
        count = VTKDATA_TETRA;
        break;

      case 2:
        count = VTKDATA_HEXAHEDRON;
        break;

      case 3:
        count = VTKDATA_WEDGE;
        break;

      case 4:
        count = VTKDATA_PYRAMID;
        break;

      default:
        fprintf(stderr, "err> %s:unknown type:%d\n", FNAME, type[ii]);
        exit(EXIT_FAILURE);
    }

#ifdef BINARY
    swap = count;
    if (isLittle) {
      BSWAP4(swap);
    }

    if (fwrite(&swap, sizeof(swap), 1, file) != 1) {
      SYS_EXIT("fwrite");
    }

#else
    if (fprintf(file, "%d\n", count) < 0) {
      SYS_EXIT("fprintf");
    }
#endif
  }

  vtk_close(file);
  free(buf);

  return;
}

/**
 * writes boundary data
 *
 * @note    none
 * @param   fileName [in] output file name
 * @param   length [in] pointer to number of faces
 * @param   array [in] pointer to node ids
 * @param   varName [in] boundary name
 * @param   num_nodes [in] pointer to total number of nodes
 * @param   isHeader [in] write header if non-zero
 * @param   len_fileName [in] length of fileName string
 * @param   len_varName [in] length of varName string
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
void
#ifdef WIN32
WRITE_VTK_BOUNDARY_BIN(const char *fileName,
                        int *length,
                        int *array,
                        const char *varName,
                        int *num_nodes,
                        int *isHeader,
                        int len_fileName,
                        int len_varName)
#else
write_vtk_boundary_bin_(const char *fileName,
                        int *length,
                        int *array,
                        const char *varName,
                        int *num_nodes,
                        int *isHeader,
                        int len_fileName,
                        int len_varName)
#endif
{
#undef FNAME
#define FNAME "write_vtk_boundary_bin_"

  FILE *file;
  char *buf;
  int ii, jj;
  
  buf = vtk_malloc(len_fileName+1);
  strncpy(buf, fileName, len_fileName);
  RSTRIP(buf);

#ifdef VERBOSE
  fprintf(stdout, "dbg> %s: opening '%s'\n", FNAME, buf);
#endif

  file = vtk_open(buf, "a+b");

  if (putc('\n', file) == EOF) {
    SYS_EXIT("putc");
  }

  if (*isHeader) {
    if (fprintf(file, "%s %d\n", VTKDATA_POINT_DATA, *num_nodes) < 0) {
      SYS_EXIT("fprintf");
    }
  }

  free(buf);
  buf = vtk_malloc(len_varName+1);
  strncpy(buf, varName, len_varName);
  RSTRIP(buf);

  fprintf(stdout, "dbg> %s: writing '%s'\n", FNAME, buf);

  if (fprintf(file, "\n%s %s char\n", VTKDATA_SCALAR, buf) < 0) {
    SYS_EXIT("fprintf");
  }

  if (fprintf(file, "%s default\n", VTKDATA_LOOKUP_TABLE) < 0) {
    SYS_EXIT("fprintf");
  }

  free(buf);
  buf = vtk_malloc(*num_nodes);
  for (ii = 0; ii < *length; ii++) {
    for (jj = 0; jj < 4; jj++) {
      int node = array[jj + ii*4];
      if (node) {
        buf[node - 1] = 1;
      }
    }
  }

#ifdef BINARY
  for (ii = 0; ii < *num_nodes; ii++) {
    if (fwrite(buf + ii, sizeof(char), 1, file) != 1) {
      SYS_EXIT("fwrite");
    }
  }

#else /* for debug */
  for (ii = 0; ii < *num_nodes; ii++) {
    if (fprintf(file, "%d ", buf[ii]) < 0) {
      SYS_EXIT("fprintf");
    }
    if ((ii + 1) % 40 == 0) {
      if (putc('\n', file) == EOF) {
        SYS_EXIT("putc");
      }
    }
  }

  if (putc('\n', file) == EOF) {
    SYS_EXIT("putc");
  }
#endif

  vtk_close(file);
  free(buf);

  return;
}

/**
 * writes scalar data
 *
 * @note    none
 * @param   fileName [in] output file name
 * @param   length [in] pointer to length of array
 * @param   array [in] pointer to data
 * @param   offset [in] offset value for array
 * @param   varName [in] name of variable
 * @param   len_fileName [in] length of fileName string
 * @param   len_varName [in] length of varName string
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
void
#ifdef WIN32
WRITE_VTK_SCALAR_BIN(const char *fileName,
                      int *length,
                      double *array,
                      int *offset,
                      const char *varName,
                      int len_fileName,
                      int len_varName)
#else
write_vtk_scalar_bin_(const char *fileName,
                      int *length,
                      double *array,
                      int *offset,
                      const char *varName,
                      int len_fileName,
                      int len_varName)
#endif
{
#undef FNAME
#define FNAME "write_vtk_scalar_bin_"

  FILE *file;
  char *buf;
  float real;
  int ii;
  int isLittle = is_little_endian();
  
  buf = vtk_malloc(len_fileName+1);
  strncpy(buf, fileName, len_fileName);
  RSTRIP(buf);

#ifdef VERBOSE
  fprintf(stdout, "dbg> %s: opening '%s'\n", FNAME, buf);
#endif

  file = vtk_open(buf, "a+b");

  if (putc('\n', file) == EOF) {
    SYS_EXIT("putc");
  }

  free(buf);
  buf = vtk_malloc(len_varName+1);
  strncpy(buf, varName, len_varName);
  RSTRIP(buf);

  fprintf(stdout, "dbg> %s: writing '%s'\n", FNAME, buf);

  if (fprintf(file, "%s %s %s\n", VTKDATA_SCALAR, buf, "float") < 0) {
    SYS_EXIT("fprintf");
  }

  if (fprintf(file, "%s %s\n", VTKDATA_LOOKUP_TABLE, "default") < 0) {
    SYS_EXIT("fprintf");
  }

#ifdef BINARY
  for (ii = *offset; ii < *length + *offset; ii++) {
    real = (float)array[ii];
    if (isLittle) {
      BSWAP4(real);
    }
    if (fwrite(&real, sizeof(real), 1, file) != 1) {
      SYS_EXIT("fwrite");
    }
  }

#else
  for (ii = *offset; ii < *length + *offset; ii++) {
    if (fprintf(file, "%E\n", array[ii]) < 0) {
      SYS_EXIT("fprintf");
    }
  }
#endif

  vtk_close(file);
  free(buf);

  return;
}

/**
 * writes vector data
 *
 * @note    none
 * @param   fileName [in] output file name
 * @param   length [in] pointer to length of array
 * @param   array [in] pointer to data
 * @param   offset [in] offset value for array
 * @param   varName [in] name of variable
 * @param   len_fileName [in] length of fileName string
 * @param   len_varName [in] length of varName string
 * @return  none
 *          call exit(EXIT_FAILURE) if error occured
 */
void
#ifdef WIN32
WRITE_VTK_VECTOR_BIN(const char *fileName,
                      int *length,
                      double *array,
                      int *offset,
                      const char *varName,
                      int len_fileName,
                      int len_varName)
#else
write_vtk_vector_bin_(const char *fileName,
                      int *length,
                      double *array,
                      int *offset,
                      const char *varName,
                      int len_fileName,
                      int len_varName)
#endif
{
#undef FNAME
#define FNAME "write_vtk_vector_bin_"

  FILE *file;
  char *buf;
  int ii, jj;
  float real;
  int isLittle = is_little_endian();
  
  buf = vtk_malloc(len_fileName+1);
  strncpy(buf, fileName, len_fileName);
  RSTRIP(buf);

#ifdef VERBOSE
  fprintf(stdout, "dbg> %s: opening '%s'\n", FNAME, buf);
#endif

  file = vtk_open(buf, "a+b");

  if (putc('\n', file) == EOF) {
    SYS_EXIT("putc");
  }

  free(buf);
  buf = vtk_malloc(len_varName+1);
  strncpy(buf, varName, len_varName);
  RSTRIP(buf);

  fprintf(stdout, "dbg> %s: writing '%s'\n", FNAME, buf);

  if (fprintf(file, "%s %s %s\n", VTKDATA_VECTOR, buf, "float") < 0) {
    SYS_EXIT("fprintf");
  }

#ifdef BINARY
  for (ii = *offset; ii < *length + *offset; ii++) {
    for (jj = 0; jj < 3; jj++) {
      real = (float)array[(*length + *offset)* jj + ii];
      if (isLittle) {
        BSWAP4(real);
      }
      if (fwrite(&real, sizeof(real), 1, file) != 1) {
        SYS_EXIT("fwrite");
      }
    }
  }

#else
  for (ii = *offset; ii < *length + *offset; ii++) {
    if (fprintf(file,
                "%E %E %E\n",
                array[(*length + *offset)*0 + ii],
                array[(*length + *offset)*1 + ii],
                array[(*length + *offset)*2 + ii]) < 0) {
      SYS_EXIT("fprintf");
    }
  }
#endif

  vtk_close(file);
  free(buf);

  return;
}
/* eof */
