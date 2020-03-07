#ifndef PARTICLEFROMTREE_H
#define PARTICLEFROMTREE_H

#include "TLorentzVector.h" 
#include "TVector3.h" 

#include "TFile.h" 
#include "TTree.h" 

#include <string>

template< typename T > class Particle
{
public:
  Particle( std::string name ) ;
  
  Particle( std::string name, TTree* tree ) ;

  void setTree( TTree* tree );

  TLorentzVector getVec() const;

  std::string name_;
  
  
private:
  
  T M_ ;
  T PX_;
  T PY_;
  T PZ_;
  T PE_;
};

template< typename T > class Vertex
{
public:
  Vertex( std::string name );
  
  Vertex( std::string name, TTree* tree ) ; 
  
  void setTree( TTree* tree );

  TVector3 getPos() const ;

  std::string name_ ; 

private:

  T X_;
  T Y_;
  T Z_;
};



template< typename T > Particle<T>::Particle( std::string name ) : 
  name_( name ) {} ;

template< typename T > Particle<T>::Particle( std::string name, TTree* tree ) : 
  name_( name )
{
  setTree( tree  );
} ;

template< typename T > void Particle<T>::setTree( TTree* tree )
{
  
  tree->SetBranchStatus( (name_ + "_TRUEP_X").c_str(), 1 );
  tree->SetBranchStatus( (name_ + "_TRUEP_Y").c_str(), 1 );
  tree->SetBranchStatus( (name_ + "_TRUEP_Z").c_str(), 1 );
  tree->SetBranchStatus( (name_ + "_TRUEP_E").c_str(), 1 );
  
  tree->SetBranchAddress( (name_ + "_TRUEP_X").c_str(), &PX_ );
  tree->SetBranchAddress( (name_ + "_TRUEP_Y").c_str(), &PY_ );
  tree->SetBranchAddress( (name_ + "_TRUEP_Z").c_str(), &PZ_ );
  tree->SetBranchAddress( (name_ + "_TRUEP_E").c_str(), &PE_ );
  
  return ;
}

template< typename T > TLorentzVector Particle<T>::getVec() const
{
  return TLorentzVector( PX_, PY_, PZ_, PE_ );
}


template< typename T > Vertex<T>::Vertex( std::string name ) : 
  name_( name ) {} ;

template< typename T > Vertex<T>::Vertex( std::string name, TTree* tree ) : 
  name_( name )
{
  setTree( tree  );
} ;

template< typename T > void Vertex<T>::setTree( TTree* tree )
{
  
  tree->SetBranchStatus( (name_ + "_TRUEORIGINVERTEX_X").c_str(), 1 );
  tree->SetBranchStatus( (name_ + "_TRUEORIGINVERTEX_Y").c_str(), 1 );
  tree->SetBranchStatus( (name_ + "_TRUEORIGINVERTEX_Z").c_str(), 1 );
  
  tree->SetBranchAddress( (name_ + "_TRUEORIGINVERTEX_X").c_str(), &X_ );
  tree->SetBranchAddress( (name_ + "_TRUEORIGINVERTEX_Y").c_str(), &Y_ );
  tree->SetBranchAddress( (name_ + "_TRUEORIGINVERTEX_Z").c_str(), &Z_ );
  
  return ;
}

template< typename T > TVector3 Vertex<T>::getPos() const
{
  return TVector3( X_, Y_, Z_ );
}

#endif 
