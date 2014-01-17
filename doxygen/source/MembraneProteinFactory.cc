// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 MembraneProteinFactory.cc
///
/// @brief 	 MembraneProteinFactory
/// @details The membrane protein factory creates a single pose from various membrane proteins
///			 loaded on the front end and initialized as membrane proteins. This single framework
///			 will then be passed off to the MembraneHub (which coordinates I/O) and sent back to the protocol it was
///			 called from (usually in pose loading)
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinFactory_cc
#define INCLUDED_core_membrane_MembraneProteinFactory_cc

// Unit Headers
#include <core/membrane/MembraneProteinFactory.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/properties/MultiChainSpanningTopology.hh>
#include <core/membrane/properties/LipidAccInfo.hh>
#include <core/membrane/properties/MultiChainLipidAccInfo.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/geometry/MembraneResidueFactory.hh>
#include <core/membrane/geometry/EmbeddingFactory.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/membrane/kinematics/MembraneFoldTree.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.MembraneProteinFactory");

/// @brief      Membrane Protein Factory
/// @details    Initializes a pose as a membrane protein

namespace core {
namespace membrane {
    
    //// Constructors /////////////////////
    
    /// @brief   Default Constructor
    /// @details Construct a Membrane Protein Factory
    ///
    /// @return [none]
    MembraneProteinFactory::MembraneProteinFactory() :
        utility::pointer::ReferenceCount(),
        prefix_file_(""),
        fullatom_(false),
        include_lips_(false)
    {
        register_options();
        init_from_cmd();
        load_required();
    }
    
    /// @brief    Default Destructor
    /// @details
    ///
    /// @note
    MembraneProteinFactory::~MembraneProteinFactory()
    {}
    
    //// Public Member Functions ///////////////////
    
    /// @brief 	 Create Membrane Protein
    /// @details Create a membrne proteins from a series of loaded membrane proteins
    ///
    /// @return  Pose (as starting structure)
    core::pose::PoseOP
    MembraneProteinFactory::create_membrane_pose() {
        
        // Create new pose, build and return
        core::pose::PoseOP pose = new core::pose::Pose();
        build_pose(pose);
        
        return pose;
    }

    //// Private Member Functions //////////////////
    
    /// @brief Register Options
    /// @details Register commandline options and flags relevant to membrane proteins
    /// these include non-major resource options
    ///
    /// @note namespaces in, membrane
    void
    MembraneProteinFactory::register_options() {
        
        using namespace basic::options;
        
        option.add_relevant( OptionKeys::in::membrane );
        option.add_relevant( OptionKeys::in::file::fullatom );
        option.add_relevant( OptionKeys::in::file::membrane_chains );
    }
    
    /// @brief      Init from CMD
    /// @details    Initialize parameters relevant to membrane options from the command
    ///             line
    ///
    /// @throws     Argument exception if chain list not specified (also this is well docuemnted, no excuse)
    void
    MembraneProteinFactory::init_from_cmd() {
        
        using namespace basic::options;
        using namespace core::membrane::util;
        
        // Ensure -in:membrane is set
        if ( option[ OptionKeys::in::membrane ].user() ) {
            if ( ! option[ OptionKeys::in::membrane ]() ) {
                throw new EXCN_NonMembrane("Cannot construct a membrane pose if the commandline option -in:membrane not passed!");
            }
        }
        
        // Initialize fullatom param
        if ( option[ OptionKeys::in::file::fullatom ].user() ) {
            fullatom_ = option[ OptionKeys::in::file::fullatom ]();
        }
        
        // Initialize prefix list from chains
        if ( option[ OptionKeys::in::file::membrane_chains ].user() ) {
            
            // Grab file and create stream
            std::string infile = option[ OptionKeys::in::file::membrane_chains ]();
            
            std::string line;
            utility::io::izstream stream (infile);
            
            std::string desc;
            
            stream.open(infile);
            if (stream) {
                
                // Grab the first line
                getline(stream, line);
                
                core::Size i = 1;
                while ( !stream.eof() ) {
                    
                    std::istringstream l(line);
                    l >> desc;
                    chains_map_.insert( std::pair< core::Size, std::string >( i, desc ) );
                    getline(stream, line);
                    i++;
                }
                
            } else {
                throw new EXCN_Illegal_Arguments("Cannot open file " + infile );
            }
        } else {
            throw new EXCN_Illegal_Arguments("Must specify a list of membrane chains");
        }
    }
    
    /// @brief Build Pose
    /// @details Create pose containing membrane/embedding residues from multi-
    /// chain input.
    /// @throws EXCN_Resource_Manager, EXCN_Membrane_Bounds
    void
    MembraneProteinFactory::build_pose( core::pose::PoseOP pose ) {
        
        using namespace core::conformation;
        using namespace core::kinematics;
        using namespace core::membrane::geometry;
        using namespace core::pose;
        
        // Keep Track of Data for the foldtree
        utility::vector1< std::pair< int, int > > peptide_edges;
        utility::vector1< std::pair< int, int > > membrane_edges;
        int nchains = (int) chains_.size();
        
        // Resize by chain size
        peptide_edges.resize(nchains);
        membrane_edges.resize(nchains);
        
        // Create multi chain pose and write teh chain map
        for ( int i = 1; i <= nchains; i++ ) {
            
            // Grab pose by chain
            core::pose::PoseOP temp_pose = chains_[i];
            
            // Add my chain to the next pose
            append_pose_to_pose(*pose, *temp_pose, true);
            
            // Add data to foldtree map
            peptide_edges[i] = std::pair< int, int >( pose->conformation().chain_begin(i), pose->conformation().chain_end(i) );
        }
        
        // Create and add Membrane residue with defaults
        core::Vector center(0, 0, 0);
        core::Vector normal(0, 0, 1);
        core::Real depth = 30.0;
        
        mrf_.add_membrane_residue(center, normal, depth, *pose, fullatom_);
        int root = (int) pose->total_residue();
        
        // For each pose/prefix pair, create an embedding residue, add it to the existing pose
        // chain and then append that chain to the master pose
        for ( core::Size i = 1; i <= 1; i++ ) {
            
            // Grab jump from initialized chain map
            core::Size jump = peptide_edges[i].first;
            
            // Grab Resources for the given chain
            core::membrane::properties::SpanningTopologyOP topology = topologies_[i];
            core::membrane::util::EmbedConfigInfoOP def = embeddings_[i];
            
            // Create a factory and apply residue
            EmbeddingFactoryOP factory = new EmbeddingFactory(pose, def, topology);
            factory->create_and_add_embedding(fullatom_, jump);
            
            // Map to membrane edges
            membrane_edges[i] = std::pair< int, int >( jump, root+i );
        }
            
        // Showing the fold tree
        TR << "Showing my current fold tree prior to integration" << std::endl;
        pose->fold_tree().show(std::cout);

        TR << "Reordering my fold tree to make the root a membrane residue" << std::endl;
        FoldTree ft( pose->fold_tree() );
        //ft.reorder( 81 );
        pose->conformation().fold_tree( ft );
        
        TR << "SHowing my new reordered fold tree" << std::endl;
        pose->fold_tree().show(std::cout);
        
        // Initialize Spanning Topology
        initialize_topology(*pose);
        
        // Should initialize multi chain lipds here
        initialize_lips_exp(*pose);
        
        // Integrate foldtree;
        integrate_foldtree( *pose, peptide_edges, membrane_edges, root );
        
        // Done!
        return;
    }
    
    /// @brief Integrate FoldTree into membrane pose
    /// @details Add required jumps between root and embedding residues
    void
    MembraneProteinFactory::integrate_foldtree(
                                               core::pose::Pose & pose,
                                               utility::vector1< std::pair< int, int > > peptide_edges,
                                               utility::vector1< std::pair< int, int > > membrane_edges,
                                               int root )
    {
        using namespace core::kinematics;
        using namespace core::membrane::kinematics;
        
        // Construct membrane foldtree from existing foldtree
        FoldTreeOP mpft = new MembraneFoldTree( peptide_edges, membrane_edges, root );
        mpft->show(std::cout);
        pose.fold_tree( *mpft );
        
    }
    
    /// @brief Initialize Spanning Topology
    /// @details Initialize spanning topology in the final pose
    void
    MembraneProteinFactory::initialize_topology( core::pose::Pose & pose )
    {
        using namespace core::membrane::properties;
        
        MultiChainSpanningTopologyOP sp = new MultiChainSpanningTopology();
        
        for ( core::Size i = 1; i <= topologies_.size(); i++ )
        {
            sp->add_by_chain( std::pair< int, SpanningTopology& >( (int) i, *topologies_[i] ) );
        }
        
        // Add to pose cache
        pose.data().set( core::pose::datacache::CacheableDataType::MULTICHAIN_SPANNING_TOPOLOGY, sp );
    }
    
    /// @brief Initialize Lipds Exposure Data
    /// @details Initialize lipid exposure data in the final pose
    void
    MembraneProteinFactory::initialize_lips_exp( core::pose::Pose & pose ) {
        
        using namespace core::membrane::properties;
        
        if ( include_lips_ ) {
            
            MultiChainLipidAccInfoOP lp = new MultiChainLipidAccInfo();
            
            for ( core::Size i = 1; i <= lipid_acc_.size(); i++ )
            {
                lp->add_by_chain( std::pair< int, LipidAccInfo& >( (int) i, *lipid_acc_[i] ) );
            }
            
            // Add to pose cache
            pose.data().set( core::pose::datacache::CacheableDataType::MULTICHAIN_LIPID_ACC_DATA, lp );
        }
    }
    
    /// @brief Load required
    /// @details Load required resources for initializing a membrane protein
    /// should make thngs thread safe (load all at once?)
    /// @note Precondition: Initialized prefix list
    ///
    /// @throws EXCN_Resource_Manager (missing reuqired resource)
    void
    MembraneProteinFactory::load_required() {
        
        using namespace basic::resource_manager;
        using namespace core::membrane::util;
        
        // Resize maps based on chains_map size
        chains_.resize(chains_map_.size());
        topologies_.resize(chains_map_.size());
        embeddings_.resize(chains_map_.size());
        lipid_acc_.resize(chains_map_.size());
        
        for ( core::Size i = 1; i <= chains_map_.size(); i++ ) {
            
            // Prefix Based Resource Tags
            std::string base_desc = chains_map_.at(i);
            std::string topo_desc = base_desc + "_span";
            std::string embed_desc = base_desc + "_embed";
            std::string lipid_desc = base_desc + "_lips";
            
            // Get pose from resource manager
            if ( ! ResourceManager::get_instance()->has_resource_with_description( base_desc ) )
            {
                throw EXCN_Resource_Definition( "Cannot load chain with description " + base_desc );
            }
            chains_[i] = basic::resource_manager::get_resource< core::pose::Pose >( base_desc );
            
            // Get topology from resource manager
            if ( ! ResourceManager::get_instance()->has_resource_with_description( topo_desc ) )
            {
                throw EXCN_Resource_Definition( "Cannot load topology with description " + topo_desc );
            }
            topologies_[i] = basic::resource_manager::get_resource< core::membrane::properties::SpanningTopology >( topo_desc );
            
            // Get embedding from resource manager
            if ( ! ResourceManager::get_instance()->has_resource_with_description( embed_desc) )
            {
                throw EXCN_Resource_Definition( "Cannot load membrane embedding with description " + embed_desc );
            }
            embeddings_[i] = basic::resource_manager::get_resource< core::membrane::util::EmbedConfigInfo >( embed_desc );
            
            // If user specified to include lips, load lipid acc data
            if ( include_lips_ ) {
                
                // Get lipids accessibility data from the resource manager
                if ( ! ResourceManager::get_instance()->has_resource_with_description( lipid_desc) )
                {
                    throw EXCN_Resource_Definition( "Cannot load lipid accessibility data with description " + lipid_desc );
                }
                lipid_acc_[i] = basic::resource_manager::get_resource< core::membrane::properties::LipidAccInfo >( lipid_desc );
            }
        }
        return;
    }

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneProteinFactory_cc


