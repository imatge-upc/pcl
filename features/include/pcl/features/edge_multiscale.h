/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2009, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#ifndef PCL_EDGE_MULTISCALE_H_
#define PCL_EDGE_MULTISCALE_H_

#include <pcl/point_types.h>
#include <pcl/features/feature.h>
#include <pcl/filters/filter.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>

namespace pcl
{
  /** \brief SHOTEstimationOMP estimates the Signature of Histograms of OrienTations (SHOT) descriptor for a given point cloud dataset
    * containing points and normals, in parallel, using the OpenMP standard.
    *
    * The suggested PointOutT is pcl::SHOT352.
    *
    * \note If you use this code in any academic work, please cite:
    *
    *   - F. Tombari, S. Salti, L. Di Stefano
    *     Unique Signatures of Histograms for Local Surface Description.
    *     In Proceedings of the 11th European Conference on Computer Vision (ECCV),
    *     Heraklion, Greece, September 5-11 2010.
    *   - F. Tombari, S. Salti, L. Di Stefano
    *     A Combined Texture-Shape Descriptor For Enhanced 3D Feature Matching.
    *     In Proceedings of the 18th International Conference on Image Processing (ICIP),
    *     Brussels, Belgium, September 11-14 2011.
    *
    * \author Samuele Salti
    * \ingroup features
    */
  template <typename PointInT, typename PointOutT = pcl::PointXYZI, typename PointInternalT = pcl::PointXYZI>
  class MultiscaleEdgeDetector;


  template <typename PointInT, typename PointOutT, typename PointInternalT>
  class MultiscaleEdgeDetector
  {
    public:
      
      /** \brief Empty constructor. */
      MultiscaleEdgeDetector (int min_K=0, int max_K=100,int max_K_color=100, int step_K=1, bool single_k=false,
                              int iters_stop_nostep=50, int iters_stop_nocolor=100,
                              float weight_level_thresh=1, float weight_level_color_thresh=1) : 
                            _min_K (min_K), _max_K (max_K), _max_K_color(max_K_color), _step_K(step_K), _single_k(single_k),
                            _iters_stop_nostep(iters_stop_nostep), _iters_stop_nocolor(iters_stop_nocolor),
                            _weight_level_thresh(weight_level_thresh), _weight_level_color_thresh(weight_level_color_thresh),
                              _kdtree_set(false)
      { };
      /** \brief Initialize the scheduler and set the number of threads to use.
        * \param nr_threads the number of hardware threads to use (0 sets the value back to automatic)
        */
      inline void
      setItersParameters (int min_K=0, int max_K=100,int max_K_color=100, int step_K=1) {
        _min_K = min_K;
        _max_K = max_K;
        _max_K_color = max_K_color;
        _step_K = step_K;
      }

      inline void
      setItersStopParameters (int iters_stop_nostep=50, int iters_stop_nocolor=100) {
        _iters_stop_nostep = iters_stop_nostep;
        _iters_stop_nocolor = iters_stop_nocolor;
      }
      inline void
      setSingleK (bool single_k) {
        _single_k = single_k;
      }

      inline void
      setThresholds (float weight_level_thresh=1, float weight_level_color_thresh=1) {
        _weight_level_thresh = weight_level_thresh;
        _weight_level_color_thresh = weight_level_color_thresh;
      }
      inline void
      setKdtree (pcl::KdTreeFLANN<PointInternalT> &kdtree) {
        _kdtree = kdtree;
        _kdtree_set = true;
      }
      void
      setInputCloud (pcl::PointCloud<PointInT> &cloud);

      void
      estimate (pcl::PointCloud<PointOutT> &output);

      float
      estimatePoint (PointInternalT &point);

    protected:
      inline void
      computeKdtree(){
         if (_cloud_set) {
         _kdtree.setInputCloud(_cloud);
         _kdtree_set = true;
         std::cout << "Kdtree computed " << std::endl;
        }
        else {

        }
      }

    private:

      int _min_K, _max_K, _max_K_color, _step_K;
      bool _single_k;
      int _iters_stop_nostep, _iters_stop_nocolor;
      float _weight_level_thresh, _weight_level_color_thresh;
      bool _kdtree_set; typename pcl::KdTreeFLANN<PointInternalT> _kdtree;
      bool _cloud_set; typename pcl::PointCloud<PointInternalT>::Ptr _cloud;
  };

}

#ifdef PCL_NO_PRECOMPILE
#include <pcl/features/impl/edge_multiscale.hpp>
#endif

#endif  //#ifndef PCL_EDGE_MULTISCALE_H_
