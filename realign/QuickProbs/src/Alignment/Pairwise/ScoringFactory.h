#pragma once
#include <string>
#include <memory>
#include <stdexcept>

#include "AminoAcidMatrices.h"
#include "NucleotideMatrices.h"

template <class T>
class ScoringFactory
{
public:
	
	static std::unique_ptr<Scoring<T>> create(std::string name, T ge, T gi) {
		auto obj = ScoringFactory<T>::create(name);
		obj->ge = ge;
		obj->gi = gi;
		return obj;
	}
	
	static std::unique_ptr<Scoring<T>> create(std::string name) {
		Scoring<T> *obj;  

		if		(name == "Blosum45")	obj = new Blosum45<T>();	
		else if (name == "Blosum50")	obj = new Blosum50<T>();
		else if (name == "Blosum62")	obj = new Blosum62<T>();
		else if (name == "Blosum80")	obj = new Blosum80<T>();
		else if (name == "Pam30")		obj = new Pam30<T>();
		else if (name == "Pam70")		obj = new Pam70<T>();
		else if (name == "Pam120")		obj = new Pam120<T>();
		else if (name == "Pam250")		obj = new Pam250<T>();
		else if (name == "Gonnet160")	obj = new Gonnet160<T>();
		else if (name == "Vtml200")		obj = new Vtml200<T>();
		else if (name == "Miqs")		obj = new Miqs<T>();
		else if (name == "MiqsFP")		obj = new MiqsFP<T>();
		else if (name == "Hoxd70")		obj = new Hoxd70<T>(); 
		else throw 
			std::runtime_error("ScoringFactory::create(): Unknown substitution matrix type.");

		return std::unique_ptr<Scoring<T>>(obj);
	}
};
