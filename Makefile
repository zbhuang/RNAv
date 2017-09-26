#CC=/opt/sfw/gcc-2/bin/gcc
CC=g++
CFLAGS=-g -c -Wall -O3
#LDLIBS=-lstdc++

all: rnav

rnav: rnav.o DynamicP.o StructureSearch.o CandScanner.o CYK.o Viterbi.o FViterbi.o ScfgBuilder.o HMMBuilder.o DataLoader.o GenomeData.o timecalc.o common.o FilterSelection.o PeakFinding.o Score.o SearchHit.o SearchOutput.o
	$(CC) rnav.o DynamicP.o StructureSearch.o CandScanner.o CYK.o Viterbi.o FViterbi.o ScfgBuilder.o HMMBuilder.o DataLoader.o GenomeData.o timecalc.o common.o FilterSelection.o PeakFinding.o Score.o SearchHit.o SearchOutput.o -o rnav

rnav.o: rnav.cpp
	$(CC) $(CFLAGS) rnav.cpp

DynamicP.o: DynamicP.cpp
	$(CC) $(CFLAGS) DynamicP.cpp 

StructureSearch.o: StructureSearch.cpp
	$(CC) $(CFLAGS) StructureSearch.cpp 

CandScanner.o: CandScanner.cpp
	$(CC) $(CFLAGS) CandScanner.cpp

CYK.o: CYK.cpp
	$(CC) $(CFLAGS) CYK.cpp

Viterbi.o: Viterbi.cpp
	$(CC) $(CFLAGS) Viterbi.cpp

FViterbi.o: FViterbi.cpp
	$(CC) $(CFLAGS) FViterbi.cpp

ScfgBuilder.o: ScfgBuilder.cpp
	$(CC) $(CFLAGS) ScfgBuilder.cpp

HMMBuilder.o: HMMBuilder.cpp
	$(CC) $(CFLAGS) HMMBuilder.cpp

DataLoader.o: DataLoader.cpp
	$(CC) $(CFLAGS) DataLoader.cpp

GenomeData.o: GenomeData.cpp
	$(CC) $(CFLAGS) GenomeData.cpp

timecalc.o: timecalc.cpp
	$(CC) $(CFLAGS) timecalc.cpp

common.o: common.cpp
	$(CC) $(CFLAGS) common.cpp

FilterSelection.o: FilterSelection.cpp
	$(CC) $(CFLAGS) FilterSelection.cpp

PeakFinding.o: PeakFinding.cpp
	$(CC) $(CFLAGS) PeakFinding.cpp

Score.o: Score.cpp
	$(CC) $(CFLAGS) Score.cpp

SearchHit.o: SearchHit.cpp
	$(CC) $(CFLAGS) SearchHit.cpp

SearchOutput.o: SearchOutput.cpp
	$(CC) $(CFLAGS) SearchOutput.cpp
