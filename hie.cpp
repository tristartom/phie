#include "hie.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>

using namespace std;
CHie::CHie(int round, string path, bool d){
	m_nCurrentRound = round;
	dirPath = path; //+ "/generated/config_hie";
	isDebug = d;
}

BOOL CHie::Create(int nParties, const vector<int>& vParams)
{
	if (vParams.size() < (unsigned)1){
		cout << "Error! This circuit needs " << 1	<< "parameters for vector length"	<< endl;
		return FALSE;
	}
	p_nCurrentRound = m_nCurrentRound;

	//Input parameters
	p_nRep = m_nRep = vParams[0]; //Specialty bits
	p_nRefBit = m_nRefCount = vParams[1]; //Reference Bits
	mode = p_nMode = vParams[2]; //Mode
	p_gateNo = gateNo = vParams[3]; //GateNo 0:auto
	p_falseCount = fpCount = vParams[4];
	//End input parameters

	isReload = false;
	m_nNumXORs = 0;

	if(isDebug) PrintInt("FalsePositive",fpCount);

	if(mode == 0){
		if(isDebug) cout << "**** No pre-processing ***** -- Round " << m_nCurrentRound << endl;
		InitParameters(nParties, vParams);
		InitGates();
		InitLayer();
		ComputeSecure_0();
		PatchParentFields();
	}else if(mode == 1){
		if(isDebug) cout << "**** Precompute Everything ****** -- Round " << m_nCurrentRound << endl;
		InitParameters(nParties, vParams);
		InitGates();
		InitLayer();
		//if(m_nCurrentRound == 0) Precompute_Everything(nParties, vParams);
		ComputeSecure_1();
		PatchParentFields();
	}else if(mode == 2){
		if(isDebug) cout << "**** Precompute more than everything ****** -- Round " << m_nCurrentRound << endl;
		InitParameters(nParties, vParams);
		Precompute_MoreThanEverything(nParties, vParams);
	}else if(mode == 3){
		if(isDebug) cout << "**** Precompute someting ****** -- Round " << m_nCurrentRound << endl;
		InitParameters(nParties, vParams);
		Precompute_Something(nParties, vParams);
	}else if(mode == 12){
		//This is a Precompute more than everything after new files are generated
		if(isDebug) cout << "**** Precompute more than everything 12 ****** -- Round " << m_nCurrentRound << endl;
		InitParameters(nParties, vParams);
		InitGates();
		InitLayer();
		ComputeSecure_2();
		//m_vOutputStart.resize(nParties, m_vInputStart[0]);
		//m_vOutputEnd.resize(nParties,m_vInputEnd[0]);

		//if(isDebug) cout << "Mode 12 Before Patch\n";
		PatchParentFields();
	}else if(mode == 13){
		//This is a Precompute something after new files are generated
		if(isDebug) cout << "**** Precompute someting 13 ****** -- Round " << m_nCurrentRound << endl;
		InitParameters(nParties, vParams);
		InitGates();
		InitLayer();
		ComputeSecure_3();

		PatchParentFields();
	}else if(mode == 100){
        if(isDebug) cout << "Merge input m1 " << endl;
        InitParameters(nParties, vParams);
        InitGates();

        vector<int> o;
        for(int i=0;i<nParties;i++){
                o.push_back( PutXORGate(m_vInputStart[i],0) );
        }

        m_vOutputStart.resize(nParties,o.front());
        m_vOutputEnd.resize(nParties,o.back());

        p_lastGateId = o.back();
	
        PatchParentFields();
  }


	else{
		cout << "Mode: " << mode << " (There is no this mode)\n";
		cout << "--------------- Crete Circuit Completed ------------------\n";
		return FALSE;
	}

	if((mode != 2) && (mode != 3) && (isDebug)){
		cout << "++++++++++++++++++++++ Summary ++++++++++++++++++++++\n";
		cout << "Esimate NonInput gates " << p_totalGates << endl;
		PrintVector("Output- Start",m_vOutputStart);
		PrintVector("OUtput- End",m_vOutputEnd);
		cout << "++++++++++++++++++++++ Crete Circuit Completed ++++++++++++++++++\n";
	}
	//return FALSE;
	return TRUE;
}

//Secure Computing
void CHie::ComputeSecure_0(){
	vector<vector<int> > v_listDistance = CalculateDistance();
	Sorting(v_listParties, v_listDistance);

	//PatchParentFields();
}

void CHie::ComputeSecure_1(){
	vector<int> t;
	int a=m_vInputStart[0];
	for(int i=0;i<m_nNumParties;i++){
		for(int j=0;j<totalInputBits;j++){
			int start = m_vInputStart[i]+j;
			a = PutXORGate(start,a);
		}

	}


	for(int i=0;i<m_nNumParties;i++){
			int start = m_vInputStart[i]+m_nCurrentRound;
			int o = PutXORGate(start,0);
			t.push_back(o);
	}


	if(isDebug) cout << "Compute Secure 1 after read input gates " << m_nNumGates << endl;

	m_vOutputStart.resize(m_nNumParties,t.front());
	m_vOutputEnd.resize(m_nNumParties,t.back());

	if(isDebug) {
		PrintVector("o1",m_vOutputStart);
		PrintVector("o2",m_vOutputEnd);
	}
	p_lastGateId = t.back();

}

void CHie::ComputeSecure_2(){
	vector< vector<int> > v_listDistance;
	vector<int> v_distance;
	//Iterate over the distance bits
	for(int i=0;i<m_nNumParties;i++){
		int startIndex = m_vInputStart[i]+m_nRefCount;
		v_distance.clear();
		//TODO: m_ndistBits
		//for(int j=0;j<(m_ndistBits+1);j++){
		for(int j=0;j<(m_ndistBits);j++){
			int index = startIndex + j;
			v_distance.push_back(index);
		}
		v_listDistance.push_back(v_distance);
	}
	//PrintVector("v_listDistance",v_listDistance);
	//Do the sorting process
	//Debug(v_listDistance);
	Sorting(v_listParties, v_listDistance);

}

void CHie::ComputeSecure_3(){
	//Find a # of input bits from the original one.
	//+1 since there are DGs + DS
	int n_bitsPerInput = m_ndistBits;
	vector<int> distance_s;

	//Ds
	for(int i=(m_nNumParties*n_bitsPerInput);i<m_nRep;i++){
		int index = m_vInputStart[0]+m_nRefCount+i;
		distance_s.push_back(index);
	}

	distance_s = ReverseBit(distance_s);

	vector< vector<int> > v_list_Dg;

	vector< vector<int> > allDistance;

	for(int i=0;i<m_nNumParties;i++){
		int startIndex = m_vInputStart[i]+m_nRefCount;
		vector<int> v_distance;
		vector<int> v_distance2;
		int countDg = 0;
		int countDgBits = 0;

		for(int j=0;j<(m_nRep-n_bitsPerInput);j++){
			int index = startIndex + j;

			if(countDgBits < n_bitsPerInput){
				v_distance.push_back(index);
				countDgBits++;
			}else{
				v_list_Dg.push_back(v_distance);
				v_distance.clear();
				v_distance.push_back(index);
				countDgBits = 1;
			}

		}
		v_list_Dg.push_back(v_distance);
		vector<int> min = FindMinValue(v_list_Dg);
		v_list_Dg.clear();
		min = Reorder(min);
		min = ReverseBit(min);

		//TODO:check output
		//Sum Ds with Dg, it function need the to be a big-endian format
		int s = PutAddGate(distance_s[0],min[0],m_ndistBits,true);
		//TODO: the lenght of vector 'sum' might be (m_ndistBits+1)
		//vector<int> sum = ReverseBit(s,m_ndistBits+1);

		int b =  m_vInputStart[i]+m_nCurrentRound;
		int result = PutMUXGate(m_nForOne2, s,b,(m_ndistBits+1),false);
		vector<int> sum = ReverseBit(result,m_ndistBits+1);

		allDistance.push_back(sum);

	}

	//PrintVector("allDistance",allDistance);
	Sorting(v_listParties, allDistance);
	//Debug(v_listParties);
	//PatchParentFields();

}

//Precompute
void CHie::Precompute_Everything(int partiesCount, vector<int> params){
	p_vPreComputeData.clear();
	int totalParties = partiesCount;
	int bits_party = GetNumBits(totalParties); //Calculate no of bits

	//Read Input from files
	//vector<int> specialty = ReadInput(totalParties,m_nRep, dirPath);

	//PrintVector("specialty",specialty);


	int total_outcomes = pow(2,partiesCount);
	int total_outcomes_bit = GetNumBits(total_outcomes);

	if(isDebug) {
		PrintInt("total_outcomes",total_outcomes);
		PrintInt("total_outcomes_bit",total_outcomes_bit);
	}
	vector<int> v_positive;
	vector<int>  v_negative;

	//Calculate the distance
	//This is a dummy data
	vector<int> t (m_nRefCount,0);
	p_vRefRead.resize(totalParties,t);
	vector<int> dist (partiesCount,0);
	vector< vector<int> > distanceAll (total_outcomes,dist);
	//if(isDebug) PrintVector("distanceAll",distanceAll);

	/*
	for(int i=0;i<total_outcomes;i++){
		//Generate I - membership vector
		vector<int> v_membership = IntToBinary(i,total_outcomes_bit);

		//Calculate distance for the given membership vector
		vector< vector<int> > v_PN = FindPositiveNegative(v_membership);
		v_positive = v_PN[0];
		v_negative = v_PN[1];

		//ComputeAllDistance Function only returns the distances of negative parties.
		//If there the size of v_positive of v_negative is zero, all elements in
		//distance vector will be -1
		vector<int> distance = ComputeAllDistance(v_positive, v_negative, specialty);

		distance = GetOutputWithFalsePositive(totalParties, distance, v_negative, 1);

		distanceAll.push_back(distance);
	}
	//End Calculate distance
	*/
	p_vPreComputeData = distanceAll;
	//Printvector("distanceAll",distanceAll);
}

void CHie::Precompute_MoreThanEverything(int partiesCount, vector<int> params){
	int totalParties = partiesCount;
	int bits_party = GetNumBits(totalParties); //Calculate no of bits

	//Read Input from files - comment if dummy
	//vector<int> specialty = ReadInput(totalParties,m_nRep, dirPath);

	//PrintVector("specialty",specialty);
	int total_outcomes = pow(2,partiesCount);
	int total_outcomes_bit = GetNumBits(total_outcomes);
	if(isDebug) {
		PrintInt("total_outcomes",total_outcomes);
		PrintInt("total_outcomes_bit",total_outcomes_bit);
	}
	vector<int> v_positive;
	vector<int>  v_negative;


	//This is a dummy data
	vector<int> t (m_nRefCount,0);
	p_vRefRead.resize(totalParties,t);
	vector<int> dist (partiesCount,0);
	vector< vector<int> > distanceAll (total_outcomes,dist);

	//if(isDebug) PrintVector("distanceAll",distanceAll);

/*
	//Calculate the distance without sorting
	vector< vector<int> > distanceAll;
	for(int i=0;i<total_outcomes;i++){
		//Generate I - membership vector
		vector<int> v_membership = IntToBinary(i,total_outcomes_bit);

		//Calculate distance for the given membership vector
		vector< vector<int> > v_PN = FindPositiveNegative(v_membership);
		v_positive = v_PN[0];
		v_negative = v_PN[1];

		//ComputeAllDistance Function only returns the distances of negative parties.
		//If there the size of v_positive of v_negative is zero, all elements in
		//distance vector will be -1
		vector<int> distance = ComputeAllDistance(v_positive, v_negative, specialty);
		distance = PutNegativeDistanceInVector(totalParties, v_negative,distance);

		distanceAll.push_back(distance);
	}
	//End Calculate distance
	*/
	//PrintVector("distanceAll",distanceAll);
	p_vPreComputeData = distanceAll;
	isReload = true;
}

void CHie::Precompute_Something(int partiesCount, vector<int> params){
	int totalParties = partiesCount;
	int bits_party = GetNumBits(totalParties); //Calculate no of bits

	//This is a dummy data.
	vector<int> t (m_nRefCount,0);
	p_vRefRead.resize(totalParties,t);
	int specialtyDist = 0;
	vector<int> dummyDist (totalParties,0);
	vector<vector<int> > geographicDist (totalParties,dummyDist);
	/*
	//Read Input from files
	vector<int> specialty = ReadInput(totalParties,m_nRep, dirPath);

	//Get ref bits, p_vRefRead values are defined after readInput is called.
	vector<int> v_tmp;
	for(int i=0;i<m_nNumParties;i++){
		int ref_bit = p_vRefRead[i][m_nCurrentRound];
		v_tmp.push_back(ref_bit);
	}

	//Separate into positive and negative set.
	vector< vector<int> > v_PN = FindPositiveNegative(v_tmp);
	vector<int>  v_positive = v_PN[0];
	vector<int> v_negative = v_PN[1];

	//Find Ds and Dg according to positive and negative set.
	int specialtyDist = Compute_Ds(specialty, v_positive,v_negative);
	vector<vector<int> > geographicDist = Compute_Dg_All(specialty, v_positive,v_negative);

	//PrintInt("specialtyDist",specialtyDist);
	//PrintVector("geographicDist",geographicDist);
	*/
	p_vPreComputeData = geographicDist;
	p_nSpecialtyDist = specialtyDist;

	if(isDebug) {
		//PrintVector("p_vPreComputeData",p_vPreComputeData);
		//PrintInt("specialtyDist",specialtyDist);
	}
	isReload = true;

}

//End Precompute

vector <int> CHie::GetOutputWithFalsePositive(int totalParties, vector<int> v, vector<int> v_negative, int noFalsePositive){
	int size = v_negative.size();
	if(size == 0) return v;

	//cout << "----- False Positive: " << noFalsePositive << " ----- \n";
	//Printvector("n",v_negative);
	//Printvector("v", v);

	if(size <= noFalsePositive) {
		//cout << " --no loop\n";
		v.clear();
		v.resize(totalParties, -1);
	}
	else{
		vector<vector <int> > t;
		for(int i=0;i<v_negative.size();i++){
			vector<int> v_m;
			v_m.push_back(v[i]); //value
			v_m.push_back(v_negative[i]); //party

			t.push_back(v_m);
		}

		//Start sorting
		for(int i=0;i<t.size();i++){
			for(int j=(t.size()-1);j>i;j--){
					if(t[j][0] < t[i][0]){
						vector<int> tmp = t[j];
						t[j] = t[i];
						t[i] = tmp;
					}
			}
		}
		//Printvector("Sort", t);
		//End sorting

		int count =0;
		v.clear();
		v.resize(totalParties,-1);
		for(int i=0;i<t.size();i++){
			int index = t[i][1];
			int value = t[i][0];
			//cout << index << " : " << value << endl;
			if(count < noFalsePositive){
				v[index] = -1;
				count++;
			}else if(count >= noFalsePositive){
				v[index] = value;
			}
		}


		//Printvector("value",v);
		//cout << "-------------------------" << endl;
	}



	return v;
}


vector<int> CHie::ComputeAllDistance(vector<int> v_positive, vector<int> v_negative, vector<int> specialty){

	//cout << "***** Start Computing All Distance ***** " << endl;
	//Printvector("v_negative",v_negative);
	//Printvector("v_positive",v_positive);

	//Caculate Specialty Distance
	int specialtyDist = Compute_Ds(specialty, v_positive,v_negative);
	vector<int> geographicDist = Compute_Dg(specialty, v_positive,v_negative);

	//PrintInt("Ds",specialtyDist);
	//Printvector("Dg", geographicDist);

	for(int i=0;i<geographicDist.size();i++){
		geographicDist[i] = geographicDist[i] + specialtyDist;
	}

	//Printvector("Distance", geographicDist);
	//cout << "***** End Computing All Distance ***** \n" << endl;
	return geographicDist;
}

vector<vector <int> > CHie::FindPositiveNegative(vector<int> v_membership){
	int size = v_membership.size();
	vector<vector <int> > ret;
	vector<int> v_positive;
	vector<int> v_negative;

	//Printvector("v_membership",v_membership);
	//Separate to negative and positive group
	for(int i=0;i<size;i++){
		if(v_membership[i] == 0) v_negative.push_back(i);
		else if(v_membership[i] == 1) v_positive.push_back(i);
	}

	ret.push_back(v_positive);
	ret.push_back(v_negative);

	return ret;
}

void CHie::GenerateNewInputFiles(vector<vector <int> > data, string path){
		path = path + "/gen.txt";
		cout <<  "Gen path " << path << endl;
		const char * c = path.c_str();
		ofstream myfile;
	  myfile.open (c);
	  myfile << "Writing this to a file.\n";
	  myfile.close();
}

vector<int> CHie::PutNegativeDistanceInVector(int totalParties,vector<int> v_negative, vector<int> v_val){
	vector<int> ret (totalParties, -1);

	for(int i=0;i<v_negative.size();i++){
		int index = v_negative[i];
		int value = v_val[i];
		ret[index] = value;
	}

	return ret;
}

vector<int> CHie::Compute_Dg(vector<int> data, vector<int> v_positive, vector<int> v_negative){
		int size = data.size();
		vector<int> ret;
		//ret.resize(size,-1);
		if((v_positive.size() == 0 ) || (v_negative.size() == 0)) {
			ret.resize(size,-1);
			return ret;
		}

		//For each party in negative group, compare them with every parties in positive group
		//Then, find the min value for each negative value;
		for(int i=0;i<v_negative.size();i++){
			vector<int> tmp;

			//This loop compares a current negative with every positive parties;
			for(int j=0;j<v_positive.size();j++){
				int dist = Compute_DistanceL2(data,v_negative[i],v_positive[j]);
				tmp.push_back(dist);
			}
			//End for-loop

			//Find mininal distance among the list
			int min = -1;
			if((v_positive.size() != 0) && (v_negative.size()!= 0)){
				min = FindMinDistance(tmp);
			}
			//ret[v_negative[i]] = min;
			ret.push_back(min);
		}
		//End for-loop

		return ret;
}

vector< vector<int> > CHie::Compute_Dg_All(vector<int> data, vector<int> v_positive, vector<int> v_negative){
		int size = data.size();
		vector<vector<int> > ret (m_nNumParties, vector<int> (m_nNumParties,-1));

		if((v_positive.size() == 0 ) || (v_negative.size() == 0)) {
			//ret.resize(size,-1);
			return ret;
		}

		//For each party in negative group, compare them with every parties in positive group
		//Then, find the min value for each negative value;
		for(int i=0;i<v_negative.size();i++){
			vector<int> tmp (m_nNumParties,-1);
			int neg_index = v_negative[i];
			//This loop compares a current negative with every positive parties;
			for(int j=0;j<v_positive.size();j++){
				int pos_index = v_positive[j];
				int dist = Compute_DistanceL2(data,neg_index,pos_index);
				tmp[pos_index] = dist;
			}
			//End for-loop

			ret[neg_index] = tmp;

		}
		//End for-loop

		return ret;
}


int CHie::Compute_DistanceL2(vector<int> data, int a,int b){
	int ret = 0;
	int value_a = data[a];
	int value_b = data[b];
	vector<int> v_a = IntToBinary(value_a,10);
	vector<int> v_b = IntToBinary(value_b,10);

	//cout << "-- Compute Dg -- " << endl;

	//Distance: c^2 = (x1-y1)^2 + (x2-y2)^2
	int distance = 0;
	for(int i=0;i<v_a.size();i++){
		int v = v_a[i] - v_b[i];
		distance += pow(v,2);
		//cout << "L2: " << v_a[i] << " - " << v_b[i] << "= " << v << endl;
	}
	ret = sqrt(distance);
	//cout << "Distance L2 = " << ret << endl;
	//cout << "-----------------" << endl;

	return ret;
}

int CHie::Compute_Ds(vector<int> specialty, vector<int> v_positive, vector<int> v_negative){
	int ret=0;
	int countN =0;
	int bitsLengtht = 5;
	//cout << "-- Compute DS -- " << endl;

	//Return -1 if there is no positive parties or negative parties.
	if((v_positive.size() == 0) || (v_negative.size() == 0)){
		return 0;
	}

	int a = specialty[v_positive[0]];
	vector<int> v_pos = IntToBinary(a,bitsLengtht);

	//Union all bits in postive set
	for(int i=1;i<v_positive.size();i++){
			int b = specialty[v_positive[i]];
			vector<int> v_b = IntToBinary(b,bitsLengtht);

			for(int index=0;index < v_pos.size();index++){
				v_pos[index] = v_pos[index] | v_b[index];
			}
	}

	a = specialty[v_negative[0]];
	vector<int> v_n = IntToBinary(a,bitsLengtht);

	//Union all bits in postive set
	for(int i=1;i<v_negative.size();i++){

			int b = specialty[v_negative[i]];
			vector<int> v_b = IntToBinary(b,bitsLengtht);

			for(int index=0;index < v_n.size();index++){
				v_n[index] = v_n[index] | v_b[index];
			}
	}


	//cout << "--- Positive ended, Find Negative vs Positve ---- \n";


	//Iterate every party in negative set
	for(int i=0;i<v_negative.size();i++){
		int t = specialty[v_negative[i]];
		vector<int> v_neg = IntToBinary(t,bitsLengtht);
		//Printvector("N", v_neg);
		for(int j=0;j<v_pos.size();j++){
			if((v_neg[j]==1) && (v_pos[j]==0)) countN++;
		}
		//cout << "Count = " << countN << endl;
	}
	//cout << "Count Total = " << countN << endl;
	//cout << "--------------------\n";
	/*for(int i=0;i<v_negative.size();i++){
		int n_value = v_negative[i];
		bool isN_in_P = false;
		//Iterate every value in P
		for(int j=0;j<v_positive.size();j++){
				int p_value = v_positive[j];
				//If value of N is in P, exit the inner loop
				if(n_value == p_value){
					isN_in_P = true;
					break;
				}
		}
		//+1 if N is not in P
		if(!isN_in_P) countN++;
	}*/
	ret = countN;
	return ret;
}

int CHie::FindMinDistance(vector<int> dist){
	int min = dist[0];
	for(int i=1;i<dist.size();i++){
		int value = dist[i];
		if((value < min)) min = value;
	}
	return min;
}

vector<int> CHie::NoMPC_Sorting(vector<int> v){
	sort(v.begin(),v.end());
	/*int size = v.size();

	for(int i=0;i<size;i++){
		for(int j=(size-1);j>i;j--){
			if(v[j] < v[i]){
				int tmp = v[j];
				v[j] = v[i];
				v[i] = tmp;
			}
		}
	}*/

	return v;
}

vector<int> CHie::MakeVectorEqualSize(int totalParties, vector<int> v, vector<int> v_negative){
	if(v.size() >= totalParties) return v;
	else{
		vector<int> ret(totalParties,-1);
		for(int i=0;i<v.size();i++){
			ret[v_negative[i]] = v[i];
		}
		return ret;
	}
}

vector<int> CHie::ReadInput(int totalCount, int inputBits, string path){
	//std::string path = "./tt_sh/config_hie/inputs";
	vector<string> input;
	vector<int> specialty_value;

	for(int i=1;i<=totalCount;i++){
		std::stringstream tmp;
		tmp << i;
		std::string read_path = path+"/inputs" + tmp.str() + ".txt";
		const char * c = read_path.c_str();
		std::ifstream myfile (c);
		std::string str;
		if(myfile.is_open()){
			cout << "Preprocess Read file: " << read_path << endl;
			while (std::getline(myfile, str))
			{
				int v = ExtractData(str,inputBits, (i-1));
				input.push_back(str);
				specialty_value.push_back(v);
			}
		}else{
			cout << "Cannot open files: " << read_path << endl;
			//i--;
		}
	}

	return specialty_value;
}

int CHie::ExtractData(string data, int inputBits, int index){
	int ret = 0;
	char *cstr = new char[data.length() + 1];
	strcpy(cstr, data.c_str());
	char * pch = strtok(cstr," ");
	vector<int> v_bits;
	vector<int> v_ref;
	//TODO: fix read file, about reference bits > 1
	int count = 0;
	while(pch != NULL){
		if(count >= m_nRefCount){
			if((count-1) == (inputBits+m_nRefCount)-1) break;
			int t = atoi(pch);
			v_bits.push_back(t);
		}else{
			int t = atoi(pch);
			v_ref.push_back(t);
			//cout << "Read: " << t << endl;
		}
		pch = strtok (NULL," ");
		count++;
	}
	p_vRefRead.push_back(v_ref);
	//PrintVector("p_vRefRead",p_vRefRead);
	//Convert specialty from binary to integer
 	ret = BinaryToInt(v_bits);
	return ret;
}

int CHie::BinaryToInt(vector<int> bitVector){
	int ret = 0;
	int size = bitVector.size();

	for(int i=0;i<size;i++){
		int bit =  bitVector[size-i-1];

		if(bit != 0){
			ret += pow(2,i);
		}
	}

	return ret;
}

vector<int> CHie::IntToBinary(int value,int length){
	vector<int> ret;

	unsigned int mask = 1 << length-1;//(sizeof(int) * 8 - 1);
	//cout << value << " - ";
	for(int i = 0; i < length; i++)
	{
		 if( (value & mask) == 0 ){
			 //cout << '0' ;
			 ret.push_back(0);
		 }
		 else{
			 ret.push_back(1);
			 //cout << '1' ;
		 }
		 mask  >>= 1;
	}
	//cout << endl ;
	return ret;
}

//************************************************************

//Computation Functions
vector< vector<int> > CHie::CalculateDistance(){
	vector<vector<int> > allOutput;

	//Compute Ds
	//TODO: Find Ds, sum all positive parties by OR gates
	//Then, every negative parties XOR with above result.
	//Then, the output AND with negative
	//Then, sum all.

	//This loop separates negative and postive parites into different arrays
	//temp: positive, temp2: negative.
	vector<int> val_pos (m_nNumParties);
	vector<int> val_neg (m_nNumParties);

	vector<int> debug;
	for(int i=0;i<m_nNumParties;i++){
		int refBits = m_vInputStart[i]+m_nCurrentRound;
		debug.push_back(refBits);
		val_pos[i] = PutMUXGate(m_vInputStart[i]+m_nRefCount,m_nForZero2,refBits,m_nRep,false);
		val_neg[i] = PutMUXGate(m_nForZero2,m_vInputStart[i]+m_nRefCount,refBits,m_nRep,false);
	}

	//Reorder gates
	vector<int> temp;
	vector<int> temp2;
	vector< vector<int> > v_listDistancePos (m_nNumParties);
	vector< vector<int> > v_listDistanceNeg (m_nNumParties);

	for(int i=0;i<m_nNumParties;i++){
		temp.resize(m_nRep);
		temp2.resize(m_nRep);
		for(int j=0;j<m_nRep;j++){
			temp[j] = (val_pos[i]+j);
			temp2[j] = (val_neg[i]+j);
		}
		v_listDistancePos[i] = temp;
		v_listDistanceNeg[i] = temp2;
	}
	////////////////////////////////////////////////////////////////////////

	//This loop union the elements in positive set and the negative set separately
	vector<int> t_pos = v_listDistancePos[0];
	vector<int> t_neg = v_listDistanceNeg[0];
	for(int i=1;i<m_nNumParties;i++){
		for(int j=0;j<m_nRep;j++){
			t_pos[j] = PutORGate(t_pos[j],v_listDistancePos[i][j]);
			t_neg[j] = PutORGate(t_neg[j],v_listDistanceNeg[i][j]);
		}
	}
	///////////////////////////////////////////////////////////////////////////

	//This 2 loops find the elements that in Negative but not in Positive
	//Set of (Neg - Pos)
	vector<int> t_np (m_nRep);
	for(int i=0;i<m_nRep;i++){
		t_np[i] = PutXORGate(t_neg[i],t_pos[i]);
	}

	for(int i=0;i<m_nRep;i++){
		t_np[i] = PutANDGate(t_neg[i],t_np[i]);
	}
	//////////////////////////////////////////////////////////////////////////
	//Sum all 1-bit
	vector<int> distance_s = SumAllBits(t_np);
	distance_s = ReverseBit(distance_s);

	//Calculate Dg and by find the hamming distance and get the smallest one.
	//Add the smallest result with Ds.
	//Loop n*n
	vector< vector<int> > ret;
  for(int i=0;i<m_nNumParties;i++)
	{
		ret.resize(m_nNumParties);
		vector<int> de (m_nNumParties,0);
    for(int j=0;j<m_nNumParties;j++){
			vector<int> v_hamm;
			v_hamm.resize(m_ndistBits,1);

			if(i!=j){
				v_hamm = PutHammDistCalculator(i,j,m_nRep);
				int refBits_1 = m_vInputStart[i]+m_nCurrentRound;
				int refBits_2 = m_vInputStart[j]+m_nCurrentRound;
				int cBit_1 = PutEQGate(refBits_1,0,1); //If refBit of i equals zero
				int cBit_2 = PutEQGate(refBits_2,1,1); //If refBit of j equals one
				int bit_1_1 = PutANDGate(cBit_1,cBit_2); //If both condition satisfy
				int o = PutMUXGate(v_hamm[0],m_nForOne,bit_1_1,m_ndistBits,false);
				for(int i=0;i<m_ndistBits;i++){
					v_hamm[i] = o+i;
				}
			}
			ret[j] = v_hamm;
    }

		vector<int> min = FindMinValue(ret);
		min = Reorder(min); //Min is little-endian
		min = ReverseBit(min); //Rearrage to big-endian

		//TODO:check output
		//Sum Ds with Dg, it function need the to be a big-endian format
		int s = PutAddGate(distance_s[0],min[0],m_ndistBits,true);
		//TODO: the lenght of vector 'sum' might be (m_ndistBits+1)
		//vector<int> sum = ReverseBit(s,m_ndistBits+1);

		int b =  m_vInputStart[i]+m_nCurrentRound;
		int result = PutMUXGate(m_nForOne2, s,b,(m_ndistBits+1),false);
		vector<int> sum = ReverseBit(result,m_ndistBits+1);
		/*for(int index=0;index<(m_ndistBits+1);index++){
			sum.push_back(s+index);
		}*/
		//sum = ReverseBit


		allOutput.push_back(sum);
  }

  return allOutput;
}

//Sort small to large
void CHie::Sorting(vector< vector<int> > party, vector< vector<int> > value){
	int vectorSize = party.size();
	//vector<int> debug;
  vector< vector<int> > v_listComp;


	value = Reorder(value);
//	PrintVector("Valuevalue",value);
	value = adjustData(value,1, value[0].size());

	for(int front_index=0;front_index<vectorSize;front_index++){
		for(int back_index=(vectorSize-1);back_index>front_index;back_index--){
		//	cout << "front " << front_index << " back " << back_index <<endl;
			vector<int> v_a = ReverseBit(value[back_index-1]);
			vector<int> v_b = ReverseBit(value[back_index]);
			vector<int> p_a = party[back_index-1];
			vector<int> p_b = party[back_index];

			int swapBit = PutNBitsComparator(v_a,v_b);


			vector< vector<int> > v_swapper = PutNBitsSwapper(swapBit,v_a,v_b);
			value[back_index] = ReverseBit(v_swapper[0]);
			value[back_index-1] = ReverseBit(v_swapper[1]);
			//vector<vector <int> >().swap(v_swapper);

			vector< vector<int> > p_swapper = PutNBitsSwapper(swapBit,p_a,p_b);
			party[back_index] = p_swapper[0];
			party[back_index-1] = p_swapper[1];
			//vector<vector <int> >().swap(p_swapper);

		}


	}
	//PrintInt("Distance value",value.size());
	//PrintInt("Distance [0] value",value[0].size());
	//Debug(value);
	//cout << "End Sorting" << endl;
	value = InsertFalsePositive(value);
	//Debug(value);
	ReorderOutput(party,value);
}

void CHie::ReorderOutput(vector< vector<int> > v){
	vector<int> newV;

	for(int i=0;i<v.size();i++){
		for(int j=0;j<v[i].size();j++){
			newV.push_back(PutXORGate(v[i][j],0));
		}
	}

	m_vOutputStart.resize(m_nNumParties, newV.front());
	m_vOutputEnd.resize(m_nNumParties, newV.back());

	p_lastGateId = newV.back();
	if(isDebug) PrintInt("Last Gate ID",p_lastGateId);
}

void CHie::ReorderOutput(vector< vector<int> > v1, vector< vector<int> > v2){
  vector<int> output;
	int outputCount = v1.size(); //- m_nRefParties;

	//PrintInt("outputCount",outputCount);

  for(int i=0;i<outputCount;i++){
    for(int j=0;j<v1[i].size();j++){
      int gate = v1[i][j];
      output.push_back(PutXORGate(gate,0));
    }
  }

  for(int i=0;i<outputCount;i++){
    for(int j=0;j<v2[i].size();j++){
      int gate = v2[i][j];
      output.push_back(PutXORGate(gate,0));
    }
  }

	//TODO:Test
  m_vOutputStart.resize(m_nNumParties, output.front());
  m_vOutputEnd.resize(m_nNumParties, output.back());

	p_lastGateId = output.back();
	if(isDebug) PrintInt("Last Gate ID",p_lastGateId);
	//vector<int>().swap(output);
  //PrintVector("m_vOutputStart",m_vOutputStart);
  //PrintVector("m_vOutputEnd",m_vOutputEnd);
}

//Insert FalsePositive in MPC stage
vector< vector<int> > CHie::InsertFalsePositive(vector< vector<int> > v){

	fpCount = 3;

	for(int i=0;i<v.size();i++){
		for(int j=0;j<v[i].size();j++){
			//if(i==0) v[i][j] = PutXORGate(0,1);
			if(i<fpCount) v[i][j] = PutXORGate(0,1);
			else v[i][j] = PutXORGate(v[i][j],0);
		}
	}

	return v;
}

vector<int> CHie::ReverseBit(vector<int> v){
	vector<int> out;

	for(int i=v.size()-1; i>=0 ;i--){
		out.push_back(PutXORGate(v[i],0));
	}

	return out;
}

vector<int> CHie::ReverseBit(int start, int length){
	vector<int> out;

	for(int i=(length-1); i>=0 ;i--){
		out.push_back(PutXORGate(start+i,0));
	}

	return out;
}

vector<int> CHie::Reorder(vector<int> v){
	vector<int> newV;

	for(int j=0;j<v.size();j++){
		newV.push_back(PutXORGate(v[j],0));
	}
	return newV;
}

vector< vector<int> > CHie::Reorder(vector< vector<int> > v){
	vector< vector<int> > newV;


	for(int i=0;i<v.size();i++){
		vector<int> tmp;
		for(int j=0;j<v[i].size();j++){
			tmp.push_back(PutXORGate(v[i][j],0));
		}
		newV.push_back(tmp);
	}

	return newV;
}

//If referValue = 0, get the value of refered parties. Else, will be all 1s
//If referValue = 1, get the value of non-referred parties. Else, will be all 1s

vector< vector<int> > CHie::adjustData(vector< vector <int> > data, int referValue,int length){

	vector< vector<int> > v_listDistance;
	vector<int> value;

	v_listDistance.resize(m_nNumParties);
	v_listRefer.resize(m_nNumParties);
	value.resize(m_nNumParties);

	vector<int> debug;
	//Change input of referred parties to all 1.
	for(int i=0;i<m_nNumParties;i++){
		//The output looks like '1 0 0 0 1' depends on the referValue
		//If referValue = 0, the value of refered bit still the same
		//If referValue = 1, invese the value of referred bit
		v_listRefer[i] = PutXORGate(m_vInputStart[i]+m_nCurrentRound,referValue);
		debug.push_back(v_listRefer[i]);
		//MUX: Control bit 1->a, 0->b
		//Control bit is v_listRefer

		value[i] = PutMUXGate(data[i][0],m_nForOne,v_listRefer[i],length,false);
		if(data[i][2] == 191) {
			cout << length << " " << value[i] << endl;
		}
	}

	//Debug(debug);

	//Reorder gates
	vector<int> temp;
	for(int i=0;i<m_nNumParties;i++){
		temp.resize(length);
		for(int j=0;j<length;j++){
			temp[j] = (value[i]+j);

		}
		//if(value[i] == 292) PrintVector("tmp",temp);
		v_listDistance[i] = temp;
	}


	return v_listDistance;
}

vector<int> CHie::FindMinValue(vector<vector<int> > listValue){

	int vectorSize = listValue.size();

	//vector<int> de;

	listValue = Reorder(listValue);
	listValue = adjustData(listValue, 0,listValue[0].size());
//	cout << "Start Find Min " << vectorSize << endl;

	//Sorting

	for(int front_index=0;front_index<vectorSize;front_index++){
		for(int back_index=(vectorSize-1);back_index>front_index;back_index--){
	//		cout << "Round " << front_index << " " << back_index << endl;
			vector<int> a = ReverseBit(listValue[back_index-1]);
			vector<int> b = ReverseBit(listValue[back_index]);

			int swapBit = PutNBitsComparator(a,b);

			//de.push_back(swapBit);

			/*if((front_index==0) && (back_index==1)){
				vector< vector<int> > v_swapper = PutNBitsSwapper(swapBit,listValue[back_index-1],listValue[back_index]);

				Debug(v_swapper[1]);
			}*/
			vector< vector<int> > v_swapper = PutNBitsSwapper(swapBit,a,b);
			listValue[back_index] = ReverseBit(v_swapper[0]);
			listValue[back_index-1] = ReverseBit(v_swapper[1]);

			//vector<int>().swap(a);
			//vector<int>().swap(b);
			//vector<vector <int> >().swap(v_swapper);

		}
	}





	return listValue[0];
}
//End computations functions

//Components
void CHie::InitParameters(int partiesCount, vector<int> params){

  m_nNumParties = partiesCount;
	m_nPartyLength = GetNumBits(m_nNumParties); // -1 in ()

	if(p_nMode == 12){
		//TODO:
		//m_ndistBits = m_nRep - 1;
		m_ndistBits = m_nRep ;
	}else if(p_nMode == 13){
		m_ndistBits = m_nRep / (m_nNumParties+1);
	}
	else{
		m_ndistBits = GetNumBits(m_nRep); //+1 out()
	}

  totalInputBits = m_nPartyLength + m_nRep + m_nRefCount;
  totalOutputBits = (m_ndistBits + m_nPartyLength) * m_nNumParties;

	//
	m_vNumVBits.clear();
	m_vInputStart.clear();
	m_vInputEnd.clear();
	m_vOutputStart.clear();
	m_vOutputEnd.clear();
	//

	m_vNumVBits.resize(m_nNumParties, 1);
	m_vInputStart.resize(m_nNumParties);
	m_vInputEnd.resize(m_nNumParties);

  m_nFrontier = 2; //Input Gate starts at 2, 0 and 1 are allocated for 0,1 bit

	p_ndistBits = m_ndistBits;
	p_nPartyLength= m_nPartyLength;

	if(isDebug) {
		PrintInt("m_nCurrentRound",m_nCurrentRound);
		PrintInt("Input",m_nRep);
		PrintInt("m_nRefCount",m_nRefCount);
		PrintInt("m_nNumParties",m_nNumParties);
		PrintInt("m_nPartyLength",m_nPartyLength);
		PrintInt("m_ndistBits",m_ndistBits);
		PrintInt("totalInputBits",totalInputBits);
		PrintInt("totalOutputBits",totalOutputBits);
	}

  //Count number of gates according to number of inputs.
	v_listPartiesStart.clear();
	v_listParties.clear();
	v_listRefParties.clear();

  v_listPartiesStart.resize(m_nNumParties);
  v_listParties.resize(m_nNumParties);
	v_listRefParties.resize(m_nNumParties);
  for (int i = 0; i < m_nNumParties; i++)
  {
    m_vInputStart[i] = m_nFrontier;
		if(mode == 100) m_nFrontier = m_nFrontier+(m_nRefCount+m_nRep);
    else m_nFrontier = m_nFrontier+(m_nRefCount+m_nRep+m_nPartyLength);
    m_vInputEnd[i] = m_nFrontier - 1;

		v_listPartiesStart[i] = (m_vInputEnd[i] - m_nPartyLength)+1; //*
    v_listRefParties[i] = m_nFrontier; //*
  }

  //Put party id gate to vector
  for(int i=0;i<(m_nNumParties);i++)
  {
    for(int j=0;j<m_nPartyLength;j++){
      int g = v_listPartiesStart[i]+j;
      v_listParties[i].push_back(g);

    }
  }
  //m_othStartGate = 1-st non-input gate
  m_othStartGate = m_nFrontier;
  nonInputGates = CalculateNonInputGates();
	m_nNumGates =  (nonInputGates + m_nFrontier);

	if(isDebug) {
		PrintVector("Input", m_vInputStart);
		PrintVector("Input End",m_vInputEnd);
	  PrintInt("nonInputGates",nonInputGates);
		PrintInt("Total Estimate Gates",m_nNumGates);
		PrintInt("Last Input Gate",m_othStartGate);
	}

	int diff = 2147483648 - m_nNumGates;
	if(diff < 0) cout << "ERROR - overflow\n";

}

vector<int> CHie::PutHammDistCalculator(int partyId1, int partyId2, int inputSize){
  vector<int> ret;
  vector<int> hammDistBits;
  int noAdder = m_ndistBits;
  int start_1 = m_vInputStart[partyId1]+m_nRefCount;
  int start_2 = m_vInputStart[partyId2]+m_nRefCount;

  //Calculate Hamming distance
  for(int i=0;i<inputSize;i++)
  {
    int t = PutXORGate(start_1 + i,start_2 + i);
    hammDistBits.push_back(t);
  }

	//Count 1-bit and return
  ret = SumAllBits(hammDistBits);

  return ret;
}

vector<int> CHie::SumAllBits(vector<int> v_bits){
  int vectorSize = v_bits.size();
  int bitNum = GetNumBits(vectorSize);

  for (int i = vectorSize; i < pow(2, bitNum); i++){
    v_bits.push_back(0);
    //cout << "[b] add" << i << endl;
  }

  int newVectorSize = pow(2, bitNum);
  int length = 1;
  while (newVectorSize > 1)
  {
    vector<int> vectorTemp;
    for (int i = 0; i < newVectorSize ; i=i+2)
    {
      int out=PutAddGate(v_bits[i], v_bits[i + 1], length, true);
      vectorTemp.push_back(out);
    }
    newVectorSize = newVectorSize / 2;
    v_bits.resize(newVectorSize);
    v_bits = vectorTemp;
    length++;
  }

  vector<int> ret;

  for(int i=(bitNum-1);i>=0;i--){
		int a = PutXORGate(v_bits[0]+i,0);
    ret.push_back(a);
  }
  return ret;
}

int CHie::PutNBitsComparator(vector<int> v_a, vector<int> v_b){
  int ret = PutGTGate(v_a[0],v_b[0],v_a.size());
  return ret;
}

vector< vector<int> > CHie::PutNBitsSwapper(int swapBit, vector<int> v_a, vector<int> v_b){
	//int n= 1;//m_ndistBits
	vector<int> v_outA;
	vector<int> v_outB;
	vector< vector<int> > ret;



	for(int i=0;i<v_a.size();i++){
		vector<int> v = PutBitSwapper(swapBit, v_a[i],v_b[i]);
		v_outA.push_back(v[0]); //a
    v_outB.push_back(v[1]); //b
	}


  ret.push_back(v_outA);
  ret.push_back(v_outB);
	ret = Reorder(ret);
  return ret;
}

//End Components

//4 AND gates, 2 ORs gates
vector<int> CHie::PutBitSwapper(int swapBit, int a, int b){
  vector<int> ret;
	int newA = PutMUXGate(a,b,swapBit, 1,false);
	int newB = PutMUXGate(b,a,swapBit, 1,false);

	ret.push_back(newA);
	ret.push_back(newB);
  //int newA = PutMUXGate(a,b,)
  /*int not_swapBit = PutNOTGate(swapBit);

  int and_1 = PutANDGate(not_swapBit, a);
  int and_2 = PutANDGate(swapBit,b);

  int or_1 = PutXORGate(and_1,and_2);


  int and_3 = PutANDGate(swapBit,a);
  int and_4 = PutANDGate(not_swapBit,b);


  int or_2 = PutXORGate(and_3,and_4);

  ret.push_back(or_1); // a'
  ret.push_back(or_2); // b'*/

  return ret;
}

void CHie::InitGates(){
	if(isDebug) cout << "Init --- ";
	m_othStartGate = m_nFrontier;

	delete [] m_pGates;

	m_pGates = new GATE[m_nNumGates];
	m_nNumXORs = 0;
	GATE* gate;

	for (int i = 0; i<m_othStartGate; i++)
	{
		gate = m_pGates + i;
		gate->p_ids = 0;
		gate->p_num = 0;
		gate->left = -1;
		gate->right = -1;
		gate->type = 0;

	}
	if(isDebug) cout << " End Init ---\n";
}

void CHie::InitLayer(){
	if(isDebug) cout << "Init Layer ---- ";
	m_nForZero = 0;
	m_nForZero2 = 0;
	m_nForOne = 0;
	m_nForOne2 = 0;

	m_nForZero = m_nFrontier;
	for(int i=0;i<m_ndistBits;i++){
		PutXORGate(0,0);
	}

	m_nForZero2 = m_nFrontier;
	for(int i=0;i<m_nRep;i++){
		PutXORGate(0,0);
	}

	m_nForOne = m_nFrontier;
	for(int i=0;i<m_ndistBits;i++){
		PutXORGate(1,0);
	}

	m_nForOne2 = m_nFrontier;
	for(int i=0;i<m_ndistBits+1;i++){
		PutXORGate(1,0);
	}

	if(isDebug) cout << " ---- End init layer \n";
}

int CHie::CalculateNonInputGates(){

	int count = 100;
	if(gateNo == 0){
		if(mode == 0){
			/*
			if(m_nRep <= 35) count = count * 2.5;
			else count = count * 1.5;
			*/
			count = (-500*m_nRep)+(239.0*m_nNumParties*m_nRep) + m_nNumParties*100;
			if((m_nNumParties >= 8 ) && (m_nNumParties < 14)) count = count * 2;
			else if((m_nNumParties >= 14 ) && (m_nNumParties <= 20)) count = count + (m_nNumParties*10000);
			else if ((m_nNumParties > 20) && (m_nNumParties < 25)) count = count + (m_nNumParties*20000);
			else if ((m_nNumParties >= 25) && (m_nNumParties <= 35))  count = 1200000 + (20000*m_nNumParties);
			//count = 6000;

		}else if(mode == 1){

			if(m_nNumParties == 4) count = 10000;
			else{
				count = (2*totalInputBits * m_nNumParties);
			}
			//((m_nRep*totalInputBits) * m_nNumParties) * 1.5;
			//count = -4500+ 1200*m_nNumParties + 0.06256*m_nRep;
			//if(m_nNumParties<10) count *= 1.2;
			//count = 10000;
			//count = 5000;//m_nNumParties * totalInputBits +10000;
			//count = -3693+ 1200*m_nNumParties + 0.06256*m_nRep;
			//if((m_nNumParties < 4) && (m_nNumParties < 7)) count = m_nNumParties * totalInputBits * 2;
			//else if((m_nNumParties >= 7) && (m_nNumParties < 11 )) count = m_nNumParties * totalInputBits * totalInputBits;
			//else count = m_nNumParties * totalInputBits * totalInputBits * 4;
			//count = m_nNumParties * ((((m_vInputEnd[0] - m_vInputStart[1]) + 1 )* 2) + totalInputBits + m_nPartyLength);
		}else if(mode == 12){
			count = -3693+ 1200*m_nNumParties + 0.06256*m_nRep;
			if((m_nRep >= 50) && (m_nRep <= 10000)) count += 500;

			if(m_nNumParties>=14) count = count + 20000;
			else count = 60000;
		}else if(mode == 13){
			count = -35000+ 9837*m_nNumParties + 0.2229*m_nRep;
			count = count * 2;
			if ((m_nNumParties >= 15) && (m_nNumParties<20)) count = count *1.5;
			else if(m_nNumParties>=20) count = count *2;
		}
	}else{
		count = gateNo;
	}

	p_totalGates = count;

  return count;
}


//1 Not = 1 Gates
int CHie::PutNOTGate(int a){
  int x = PutXORGate(a,1);
  return x;
}

//1 OR = 3 Gates
int CHie::PutORGate(int a, int b){
  int x = PutXORGate(a,b);
  int y = PutANDGate(a,b);
  int z = PutXORGate(x,y);

  return z;
}

int CHie::PutNANDGate(int a, int b){
  int and_1 = PutANDGate(a,b);
  int ret = PutNOTGate(and_1);

  return ret;
}

int CHie::PutNAND2Gate(int a, int b, int c){
  int and_1 = PutANDGate(a,b);
  int and_2 = PutANDGate(and_1,c);
  int ret = PutNOTGate(and_2);

  return ret;
}

int CHie::PutNORGate(int a, int b){
  int or_1 = PutORGate(a,b);
  int ret = PutNOTGate(or_1);

  return ret;
}


int CHie::GetNumBits(int decimal){
    int num_bits = ceil(log2(decimal));
    if(decimal == 2) num_bits = 2;
		//else if(decimal == 4) num_bits = 3;

    return num_bits;
}

void CHie::PrintInt(const string& name, int i){
  cout << "|debug| " << name << " -- "<< i << endl;
}

void CHie::PrintVector(const string& name, const vector<int>& params){
  cout << "|vector| " << name << "(" << params.size() << ") :";
  for (vector<int>::const_iterator i = params.begin(); i != params.end(); ++i){
		cout << *i << " ";
	}
  cout << endl;
}

void CHie::PrintVector(const string& name, const vector< vector<int> >& v){
	cout << "|vector2| " << name << "(" << v.size() << ") " << endl;
	for(int i=0;i<v.size();i++){
		cout << "   " << i << ") " ;
		for(int j=0;j<v[i].size();j++){
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
}

void CHie::Debug(vector<int> v){
  cout << "----- Debugging -----" << endl;
  vector<int> t;
  for(int i=0;i<v.size();i++){
    int a = PutXORGate(v[i],0);
    t.push_back(a);
  }

  m_vOutputStart.resize(m_nNumParties,t.front());
  m_vOutputEnd.resize(m_nNumParties,t.back());

  PrintVector("m_vOutputStart",m_vOutputStart);
  PrintVector("m_vOutputEnd",m_vOutputEnd);


}

void CHie::Debug(vector<vector<int> > v){
	vector<int> out;
	//PrintVector("v",v[0]);
	//PrintVector("v",v[1]);
	//PrintVector("v",v[2]);
	//PrintVector("v",v[3]);
	int zero = PutXORGate(1,1);
	for(int i=0;i<v.size();i++){
		for(int j=0;j<v[i].size();j++){
			int a = PutXORGate(v[i][j],zero);
			out.push_back(a);
		}
	}

	m_vOutputStart.resize(m_nNumParties,out.front());
	m_vOutputEnd.resize(m_nNumParties,out.back());
	//m_vOutputStart.resize(m_nNumParties,v[0].front());
	//am_vOutputEnd.resize(m_nNumParties,v[0].back());


	//PrintVector("m_vOutputStart",m_vOutputStart);
	//PrintVector("m_vOutputEnd",m_vOutputEnd);

}

/*void CHie::PatchParentFields()
{
	cout << "--- Patch start -- ";

	vector<int> counter_map(m_nNumGates,0);
	GATE* gate;

	GATE* c;
	int child;

	for(int i=m_othStartGate; i<p_lastGateId; i++)
	{
		cout << i << " " << m_nNumGates << endl;
		gate = m_pGates + i;

		child = gate->left;

		if( child >= 2 )
		{
			c = m_pGates + child;
			if( counter_map[child] == 0 )
			{
				c->p_ids = New(c->p_num);
			}
			c->p_ids[ counter_map[child]++ ] = i;
		}

		child = gate->right;

		if( child >= 2 )
		{
			c = m_pGates + child;
			if( counter_map[child] == 0 )
			{
				c->p_ids = New(c->p_num);
			}
			c->p_ids[ counter_map[child]++ ] = i;
		}
	}

	cout << " Patch -- end ------ \n";
}*/
