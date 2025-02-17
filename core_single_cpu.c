/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- */

/*
 * This code has been contributed by the DARPA HPCS program.  Contact
 * David Koester <dkoester@mitre.org> or Bob Lucas <rflucas@isi.edu>
 * if you have questions.
 *
 * GUPS (Giga UPdates per Second) is a measurement that profiles the memory
 * architecture of a system and is a measure of performance similar to MFLOPS.
 * The HPCS HPCchallenge RandomAccess benchmark is intended to exercise the
 * GUPS capability of a system, much like the LINPACK benchmark is intended to
 * exercise the MFLOPS capability of a computer.  In each case, we would
 * expect these benchmarks to achieve close to the "peak" capability of the
 * memory system. The extent of the similarities between RandomAccess and
 * LINPACK are limited to both benchmarks attempting to calculate a peak system
 * capability.
 *
 * GUPS is calculated by identifying the number of memory locations that can be
 * randomly updated in one second, divided by 1 billion (1e9). The term "randomly"
 * means that there is little relationship between one address to be updated and
 * the next, except that they occur in the space of one half the total system
 * memory.  An update is a read-modify-write operation on a table of 64-bit words.
 * An address is generated, the value at that address read from memory, modified
 * by an integer operation (add, and, or, xor) with a literal value, and that
 * new value is written back to memory.
 *
 * We are interested in knowing the GUPS performance of both entire systems and
 * system subcomponents --- e.g., the GUPS rating of a distributed memory
 * multiprocessor the GUPS rating of an SMP node, and the GUPS rating of a
 * single processor.  While there is typically a scaling of FLOPS with processor
 * count, a similar phenomenon may not always occur for GUPS.
 *
 * For additional information on the GUPS metric, the HPCchallenge RandomAccess
 * Benchmark,and the rules to run RandomAccess or modify it to optimize
 * performance -- see http://icl.cs.utk.edu/hpcc/
 *
 */

/*
 * This file contains the computational core of the single cpu version
 * of GUPS.  The inner loop should easily be vectorized by compilers
 * with such support.
 *
 * This core is used by both the single_cpu and star_single_cpu tests.
 */

/* Number of updates to table (suggested: 4x number of table entries) */
#define NUPDATE (4 * TableSize)

#define POLY 0x0000000000000007UL
#define PERIOD 1317624576693539401L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>
#include <assert.h>

static double
timestamp() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (1e-6*((double)tv.tv_usec));
}

/* Utility routine to start random number generator at Nth step */
static uint64_t
HPCC_starts(int64_t n)
{
  int i, j;
  uint64_t m2[64];
  uint64_t temp, ran;

  while (n < 0) n += PERIOD;
  while (n > PERIOD) n -= PERIOD;
  if (n == 0) return 0x1;

  temp = 0x1;
  for (i=0; i<64; i++) {
    m2[i] = temp;
    temp = (temp << 1) ^ ((int64_t) temp < 0 ? POLY : 0);
    temp = (temp << 1) ^ ((int64_t) temp < 0 ? POLY : 0);
  }

  for (i=62; i>=0; i--)
    if ((n >> i) & 1)
      break;

  ran = 0x2;
  while (i > 0) {
    temp = 0;
    for (j=0; j<64; j++)
      if ((ran >> j) & 1)
        temp ^= m2[j];
    ran = temp;
    i -= 1;
    if ((n >> i) & 1)
      ran = (ran << 1) ^ ((int64_t) ran < 0 ? POLY : 0);
  }

  return ran;
}


static void
RandomAccessUpdate(uint64_t TableSize, uint64_t *Table) {
  uint64_t i;
  uint64_t ran[128];              /* Current random numbers */
  int j;

  /* Perform updates to main table.  The scalar equivalent is:
   *
   *     uint64_t ran;
   *     ran = 1;
   *     for (i=0; i<NUPDATE; i++) {
   *       ran = (ran << 1) ^ (((s64Int) ran < 0) ? POLY : 0);
   *       table[ran & (TableSize-1)] ^= ran;
   *     }
   */
  for (j=0; j<128; j++) {
    ran[j] = HPCC_starts ((NUPDATE/128) * j);
  }
  
  for (i=0; i<NUPDATE/128; i++) {
    for (j=0; j<128; j++) {
      ran[j] = (ran[j] << 1) ^ ((int64_t) ran[j] < 0 ? POLY : 0);
      Table[ran[j] & (TableSize-1)] ^= ran[j];
    }
  }
}

int main(int argc, char *argv[]) {
  uint64_t i;
  uint64_t temp;
  double realtime;              /* Real time to update table */
  double totalMem;
  uint64_t *Table;
  uint64_t logTableSize, TableSize;
  FILE *outFile = NULL;

  /* calculate local memory per node for the update table */
  totalMem = 1UL<<32;
  totalMem /= sizeof(uint64_t);

  /* calculate the size of update array (must be a power of 2) */
  for (totalMem *= 0.5, logTableSize = 0, TableSize = 1;
       totalMem >= 1.0;
       totalMem *= 0.5, logTableSize++, TableSize <<= 1)
    ; /* EMPTY */

  Table = (uint64_t*)malloc(sizeof(uint64_t)*TableSize);
  assert(Table != NULL);

  /* Print parameters for run */

  printf("Main table size   = 2 ^ %lu = %lu words\n", logTableSize,TableSize);
  printf("Number of updates = %lu\n", NUPDATE);

  /* Initialize main table */
  for (i=0; i<TableSize; i++) {
    Table[i] = i;
 }

  /* Begin timing here */
  realtime = timestamp();

  RandomAccessUpdate( TableSize, Table );

  /* End timed section */
  realtime = timestamp() - realtime;

  /* make sure no division by zero */
  double GUPs = (realtime > 0.0 ? 1.0 / realtime : -1.0);
  GUPs *= 1e-9*NUPDATE;
  /* Print timing results */
  printf( "Real time used = %.6f seconds\n", realtime);
  printf( "%.9f Billion(10^9) Updates    per second [GUP/s]\n", GUPs );
  
  /* Verification of results (in serial or "safe" mode; optional) */
  temp = 0x1;
  for (i=0; i<NUPDATE; i++) {
    temp = (temp << 1) ^ (((int64_t) temp < 0) ? POLY : 0);
    Table[temp & (TableSize-1)] ^= temp;
  }

  temp = 0;
  for (i=0; i<TableSize; i++)
    if (Table[i] != i)
      temp++;

  printf("Found %lu errors in %lu locations (%s).\n",
         temp, TableSize, (temp <= 0.01*TableSize) ? "passed" : "failed");
  

  free( Table );
  
  return 0;
}
