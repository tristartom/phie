#ifndef HIE_H_BY_BEBE_
#define HIE_H_BY_BEBE_

#include "circuit.h"
#include <cassert>

class CHie : public CCircuit
{
public:
	CHie(int round, string path, bool d);
	BOOL Create(int nParties, const vector<int>& vParams);

public:
	//Precompute
	void Precompute_Everything(int partiesCount, vector<int> params);
	void Precompute_MoreThanEverything(int partiesCount, vector<int> params);
	void Precompute_Something(int partiesCount, vector<int> params);
	vector<int> ComputeAllDistance(vector<int> v_positive, vector<int> v_negative, vector<int> specialty);


	vector<int> Compute_Dg(vector<int> data, vector<int> v_positive, vector<int> v_negative);
	vector< vector<int> > Compute_Dg_All(vector<int> data, vector<int> v_positive, vector<int> v_negative);
	int Compute_Ds(vector<int> specialty, vector<int> v_positive, vector<int> v_negative);
	int Compute_DistanceL2(vector<int> data, int a,int b);
	int FindMinDistance(vector<int> dist);
	vector<int> NoMPC_Sorting(vector<int> v);
	vector<int> MakeVectorEqualSize(int totalParties, vector<int> v, vector<int> v_negative);

	vector<int> ReadInput(int totalCount, int inputBits, string path);
	int ExtractData(string data, int inputBits, int index);
	int BinaryToInt(vector<int> bitVector);
 	vector<int> IntToBinary(int value, int length);

	vector<int> GetOutputWithFalsePositive(int totalParties, vector<int> v, vector<int> v_negative, int noFalsePositive);
	vector<vector <int> > FindPositiveNegative(vector<int> v_membership);
	vector<int> PutNegativeDistanceInVector(int totalParties,vector<int> v_negative, vector<int> v_val);

	void GenerateNewInputFiles(vector<vector <int> > data, string path);
	// End Precompute

  void InitGates();
  void InitLayer();
  void InitParameters(int partiesCount, vector<int> params);
	int GetNumBits(int);
  int CalculateNonInputGates();

  //Computations
  vector< vector<int> > CalculateDistance();
  void Sorting(vector< vector<int> > party,vector< vector<int> > value);
  void ReorderOutput(vector< vector<int> > v);
  void ReorderOutput(vector< vector<int> > v1,vector< vector<int> > v2);
  vector< vector<int> > adjustData(vector< vector <int> > data, int referValue, int length);
  vector<int> FindMinValue(vector<vector<int> > listValue);
  vector<int> ReverseBit(vector<int> v);
	vector<int> ReverseBit(int start, int length);
  vector<int> Reorder(vector<int> v);
  vector< vector<int> > Reorder(vector< vector<int> > v);
	vector< vector<int> > InsertFalsePositive(vector< vector<int> > v);
	void ComputeSecure_0();
	void ComputeSecure_1();
	void ComputeSecure_2();
	void ComputeSecure_3();

  //Components
  vector<int> PutHammDistCalculator(int partyId1, int partyId2, int inputSize);
  vector<int> SumAllBits(vector<int> v_bits);
  int PutNBitsComparator(vector<int> v_a, vector<int> v_b);
  vector<int> PutBitSwapper(int swapBit, int a, int b);
  vector< vector<int> > PutNBitsSwapper(int swapBit, vector<int> v_a, vector<int> v_b);

  //Additional gates;
  int PutNOTGate(int a);
  int PutORGate(int a, int b);
  int PutNANDGate(int a, int b);
  int PutNAND2Gate(int a, int b, int c);
  int PutNORGate(int a, int b);

  //For debugging
  void PrintInt(const string& name, int i);
  void PrintVector(const string& name, const vector<int>& params);
	void PrintVector(const string& name, const vector< vector<int> >& v);
  void Debug(vector<int> v);
  void Debug(vector<vector<int> > v);

private:
	int fpCount; //no of false positive

  int m_nRep; //number of bits represents input value
	int m_nRefCount; //number of referenced bits
  int m_nPartyLength; //number of bits represents party id
  int m_ndistBits; //number of bits represents hamming distance
  //int m_nRefParties; //number of reference parties
  int m_nSortedParties; //number of non-reference parties;
  int m_nForZero;
	int m_nForZero2;
  int m_nForOne;
	int m_nForOne2;
	int mode;

	int m_nCurrentRound;
	string dirPath;
	int gateNo;

	bool isDebug;

  int nonInputGates;
  int totalInputBits; //number of bits represent party id + value
  int totalOutputBits;
  int hammCalculatorCount; //number of hamming calculator
  int counterCount; //number of counter
  int comparatorCount; //number of comparator
  int swapperCount; //number of swapper

  vector<int> v_listRefParties;
  vector<int> v_listSortedParties;
  vector<int> v_listPartiesStart;
  vector<int> v_listRefer;
  vector< vector<int> > v_listParties;


};

#endif //HIE_H_BY_BEBE_
