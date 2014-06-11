#ifndef _FWPFCLUSTERRPZPROXYBUILDER_H_
#define _FWPFCLUSTERRPZPROXYBUILDER_H_

// -*- C++ -*-
//
// Package:     ParticleFlow
// Class  :     FWPFClusterRPZProxyBuilder, FWPFEcalClusterRPZProxyBuilder, FWPFHcalClusterRPZProxyBuilder
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Simon Harris
//

// User include files
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "Fireworks/Core/interface/FWSimpleProxyBuilderTemplate.h"
#include "Fireworks/ParticleFlow/interface/FWPFGeom.h"
#include "Fireworks/ParticleFlow/interface/FWPFClusterRPZUtils.h"

//-----------------------------------------------------------------------------
// FWPFClusterRPZProxyBuilder
//-----------------------------------------------------------------------------

class FWPFClusterRPZProxyBuilder : public FWSimpleProxyBuilderTemplate<reco::PFCluster>
{
   public:
   // ---------------- Constructor(s)/Destructor ----------------------
      FWPFClusterRPZProxyBuilder();
      virtual ~FWPFClusterRPZProxyBuilder();

   // --------------------- Member Functions --------------------------
      using FWSimpleProxyBuilderTemplate<reco::PFCluster>::build;
      virtual void build( const reco::PFCluster &iData, unsigned int iIndex, TEveElement &oItemHolder, const FWViewContext *vc );
      using FWSimpleProxyBuilderTemplate<reco::PFCluster>::scaleProduct;
      virtual void scaleProduct( TEveElementList *parent, FWViewType::EType, const FWViewContext *vc );
      using FWSimpleProxyBuilderTemplate<reco::PFCluster>::havePerViewProduct;
      virtual bool havePerViewProduct( FWViewType::EType ) const { return true; }
      using FWSimpleProxyBuilderTemplate<reco::PFCluster>::cleanLocal;
      virtual void cleanLocal() { m_clusters.clear(); }

      REGISTER_PROXYBUILDER_METHODS();

   protected:
   // ----------------------- Data Members ----------------------------
      std::vector<ScalableLines> m_clusters;
      FWPFClusterRPZUtils        *m_clusterUtils;

   // --------------------- Member Functions --------------------------
      virtual void sharedBuild( const reco::PFCluster &cluster, unsigned int iIndex, TEveElement &oItemHolder, 
                                const FWViewContext *vc, float radius );

   private:
      FWPFClusterRPZProxyBuilder( const FWPFClusterRPZProxyBuilder& );                    // Disable default
      const FWPFClusterRPZProxyBuilder& operator=( const FWPFClusterRPZProxyBuilder& );   // Disable default
};
//=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_


//-----------------------------------------------------------------------------
// FWPFEcalClusterRPZProxyBuilder
//-----------------------------------------------------------------------------

class FWPFEcalClusterRPZProxyBuilder : public FWPFClusterRPZProxyBuilder
{
   public:
   // ---------------- Constructor(s)/Destructor ----------------------
      FWPFEcalClusterRPZProxyBuilder(){}
      virtual ~FWPFEcalClusterRPZProxyBuilder(){}

   // --------------------- Member Functions --------------------------
      using FWSimpleProxyBuilderTemplate<reco::PFCluster>::build;
      virtual void build( const reco::PFCluster &iData, unsigned int iIndex, TEveElement &oItemHolder, const FWViewContext *vc );

      REGISTER_PROXYBUILDER_METHODS();

   private:
      FWPFEcalClusterRPZProxyBuilder( const FWPFEcalClusterRPZProxyBuilder& );
      const FWPFEcalClusterRPZProxyBuilder& operator=( const FWPFEcalClusterRPZProxyBuilder& );
};
//=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_


//-----------------------------------------------------------------------------
// FWPFHcalClusterRPZProxyBuilder
//-----------------------------------------------------------------------------

class FWPFHcalClusterRPZProxyBuilder : public FWPFClusterRPZProxyBuilder
{
   public:
   // ---------------- Constructor(s)/Destructor ----------------------
      FWPFHcalClusterRPZProxyBuilder(){}
      virtual ~FWPFHcalClusterRPZProxyBuilder(){}

   // --------------------- Member Functions --------------------------
      using FWSimpleProxyBuilderTemplate<reco::PFCluster>::build;
      virtual void build( const reco::PFCluster &iData, unsigned int iIndex, TEveElement &oItemHolder, const FWViewContext *vc );

      REGISTER_PROXYBUILDER_METHODS();

   private:
      FWPFHcalClusterRPZProxyBuilder( const FWPFHcalClusterRPZProxyBuilder& );
      const FWPFHcalClusterRPZProxyBuilder& operator=( const FWPFHcalClusterRPZProxyBuilder& );
};
#endif
//=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_=_
