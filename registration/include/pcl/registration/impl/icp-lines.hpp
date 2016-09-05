/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc
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
 * $Id$
 *
 */

#ifndef PCL_REGISTRATION_IMPL_ICP_LINES_HPP_
#define PCL_REGISTRATION_IMPL_ICP_LINES_HPP_

#include <pcl/registration/boost.h>
#include <pcl/correspondence.h>

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> void
pcl::IterativeClosestPointLine<PointSource, PointTarget, Scalar>::transformCloud (
    const PointCloudSource &input, 
    PointCloudSource &output, 
    const Matrix4 &transform)
{
  Eigen::Vector4f pt (0.0f, 0.0f, 0.0f, 1.0f), pt_t;
  Eigen::Matrix4f tr = transform.template cast<float> ();

  for (size_t i = 0; i < input.size (); ++i)
  {
    const uint8_t* data_in = reinterpret_cast<const uint8_t*> (&input[i]);
    uint8_t* data_out = reinterpret_cast<uint8_t*> (&output[i]);
    memcpy (&pt[0], data_in + x_idx_offset_, sizeof (float));
    memcpy (&pt[1], data_in + y_idx_offset_, sizeof (float));
    memcpy (&pt[2], data_in + z_idx_offset_, sizeof (float));

    if (!pcl_isfinite (pt[0]) || !pcl_isfinite (pt[1]) || !pcl_isfinite (pt[2])) 
      continue;

    pt_t = tr * pt;

    memcpy (data_out + x_idx_offset_, &pt_t[0], sizeof (float));
    memcpy (data_out + y_idx_offset_, &pt_t[1], sizeof (float));
    memcpy (data_out + z_idx_offset_, &pt_t[2], sizeof (float));
  }
  
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> void
pcl::IterativeClosestPointLine<PointSource, PointTarget, Scalar>::computeTransformation (
    PointCloudSource &output, const Matrix4 &guess)
{
  // Point cloud containing the correspondences of each point in <input, indices>
  PointCloudSourcePtr target_transformed (new PointCloudSource);

  nr_iterations_ = 0;
  converged_ = false;

  // Initialise final transformation to the guessed one
  final_transformation_ = guess.inverse().eval();

  // If the guessed transformation is non identity
  if (guess != Matrix4::Identity ())
  {
    target_transformed->resize (target_->size ());
     // Apply guessed transformation prior to search for neighbours
    transformCloud (*target_, *target_transformed, final_transformation_);
  }
  else
    copyPointCloud (*target_, *target_transformed); 
 
  transformation_ = Matrix4::Identity ();

  correspondence_estimation_->setInputTarget (input_);

  convergence_criteria_->setMaximumIterations (max_iterations_);
  convergence_criteria_->setRelativeMSE (euclidean_fitness_epsilon_);
  convergence_criteria_->setTranslationThreshold (transformation_epsilon_);
  convergence_criteria_->setRotationThreshold (1.0 - transformation_epsilon_);
  PointCloudSourcePtr target_shrinked (new PointCloudSource);

  // Repeat until convergence
  int idx = 0;
  do
  {
    idx++;
    previous_transformation_ = transformation_;
    std::cout << "1" << std::endl;
    shrinkPointCloud(*target_transformed, *target_shrinked, input_->points[0].z);
    correspondence_estimation_->setInputSource (target_transformed);
    // correspondence_estimation_->setInputSource (target_shrinked);
    pcl::io::savePLYFileASCII ("target_shrinked_"+boost::lexical_cast<std::string>(idx)+".ply", *target_shrinked);
    pcl::io::savePLYFileASCII ("target_transformed_"+boost::lexical_cast<std::string>(idx)+".ply", *target_transformed);
    pcl::io::savePLYFileASCII ("input"+boost::lexical_cast<std::string>(idx)+".ply", *input_);
    std::cout << "2" << std::endl;

    if (use_reciprocal_correspondence_)
      correspondence_estimation_->determineReciprocalCorrespondences (*correspondences_, corr_dist_threshold_);
    else
      correspondence_estimation_->determineCorrespondences (*correspondences_, corr_dist_threshold_);
    std::cout << "3" << std::endl;

    for (int i = 0; i < (*correspondences_).size(); ++i)
    {
      // PCL_ERROR("Point correspondence: Source %lu target %lu \n", 
      //   (*correspondences_)[i].index_query, 
      //   (*correspondences_)[i].index_match);
    }
    // if (correspondence_rejectors_.empty ())
    // ToDO: Correspondece rejectors?!!
    // CorrespondencesPtr temp_correspondences (new Correspondences (*correspondences_));
    // for (size_t i = 0; i < correspondence_rejectors_.size (); ++i)
    // {
    //   PCL_ERROR("There is a correspondence rejector");
    //   registration::CorrespondenceRejector::Ptr& rej = correspondence_rejectors_[i];
    //   PCL_DEBUG ("Applying a correspondence rejector method: %s.\n", rej->getClassName ().c_str ());
    //   // if (rej->requiresSourcePoints ())
    //   //   rej->setSourcePoints (input_transformed_blob);
    //   // if (rej->requiresSourceNormals () && source_has_normals_)
    //   //   rej->setSourceNormals (input_transformed_blob);
    //   rej->setInputCorrespondences (temp_correspondences);
    //   rej->getCorrespondences (*correspondences_);
    //   // Modify input for the next iteration
    //   if (i < correspondence_rejectors_.size () - 1)
    //     *temp_correspondences = *correspondences_;
    // }

    size_t cnt = correspondences_->size ();
    // Check whether we have enough correspondences
    if (static_cast<int> (cnt) < min_number_correspondences_)
    {
      PCL_ERROR ("[pcl::%s::computeTransformation] Not enough correspondences found. Relax your threshold parameters.\n", getClassName ().c_str ());
      convergence_criteria_->setConvergenceState(pcl::registration::DefaultConvergenceCriteria<Scalar>::CONVERGENCE_CRITERIA_NO_CORRESPONDENCES);
      converged_ = false;
      break;
    }

    // Estimate the transform
    transformation_estimation_->estimateRigidTransformation (*target_transformed, *input_, *correspondences_, transformation_);
    // Tranform the data
    transformCloud (*target_transformed, *target_transformed, transformation_);

    // Obtain the final transformation    
    final_transformation_ = transformation_ * final_transformation_;
    std::cout << transformation_ << std::endl;
    ++nr_iterations_;
    std::cout << "h" << std::endl;
    converged_ = static_cast<bool> ((*convergence_criteria_));
    std::cout << "he" << std::endl;

  }
  while (!converged_);
  final_transformation_ = final_transformation_.inverse().eval();

  // Transform the input cloud using the final transformation
  PCL_DEBUG ("Transformation is:\n\t%5f\t%5f\t%5f\t%5f\n\t%5f\t%5f\t%5f\t%5f\n\t%5f\t%5f\t%5f\t%5f\n\t%5f\t%5f\t%5f\t%5f\n", 
      final_transformation_ (0, 0), final_transformation_ (0, 1), final_transformation_ (0, 2), final_transformation_ (0, 3),
      final_transformation_ (1, 0), final_transformation_ (1, 1), final_transformation_ (1, 2), final_transformation_ (1, 3),
      final_transformation_ (2, 0), final_transformation_ (2, 1), final_transformation_ (2, 2), final_transformation_ (2, 3),
      final_transformation_ (3, 0), final_transformation_ (3, 1), final_transformation_ (3, 2), final_transformation_ (3, 3));

  transformCloud (*input_, output, final_transformation_);
  // shrinkPointCloud(*target_, output, input_->points[0].z);
}

template <typename PointSource, typename PointTarget, typename Scalar> void
pcl::IterativeClosestPointLine<PointSource, PointTarget, Scalar>::determineRequiredBlobData ()
{
  need_source_blob_ = false;
  need_target_blob_ = false;
  // Check estimator
  need_source_blob_ |= correspondence_estimation_->requiresSourceNormals ();
  need_target_blob_ |= correspondence_estimation_->requiresTargetNormals ();
  // Add warnings if necessary
  if (correspondence_estimation_->requiresSourceNormals ())
  {
      PCL_WARN("[pcl::%s::determineRequiredBlobData] Estimator expects source normals, but we can't provide them.\n", getClassName ().c_str ());
  }
  if (correspondence_estimation_->requiresTargetNormals () && !target_has_normals_)
  {
      PCL_WARN("[pcl::%s::determineRequiredBlobData] Estimator expects target normals, but we can't provide them.\n", getClassName ().c_str ());
  }
  // Check rejectors
  for (size_t i = 0; i < correspondence_rejectors_.size (); i++)
  {
    registration::CorrespondenceRejector::Ptr& rej = correspondence_rejectors_[i];
    need_source_blob_ |= rej->requiresSourcePoints ();
    need_source_blob_ |= rej->requiresSourceNormals ();
    need_target_blob_ |= rej->requiresTargetPoints ();
    need_target_blob_ |= rej->requiresTargetNormals ();
    if (rej->requiresSourceNormals () && !source_has_normals_)
    {
      PCL_WARN("[pcl::%s::determineRequiredBlobData] Rejector %s expects source normals, but we can't provide them.\n", getClassName ().c_str (), rej->getClassName ().c_str ());
    }
    if (rej->requiresTargetNormals () && !target_has_normals_)
    {
      PCL_WARN("[pcl::%s::determineRequiredBlobData] Rejector %s expects target normals, but we can't provide them.\n", getClassName ().c_str (), rej->getClassName ().c_str ());
    }
  }
}
      
template <typename PointSource, typename PointTarget, typename Scalar> void
pcl::IterativeClosestPointLine<PointSource, PointTarget, Scalar>::shrinkPointCloud (
    const PointCloudSource &points, PointCloudSource & points_shrinked, Scalar depth)
{
    points_shrinked.clear();
  for (size_t i = 0; i < points.points.size (); ++i){
        PointSource p;
      p.x = points.points[i].x*depth/points.points[i].z;
      p.y = points.points[i].y*depth/points.points[i].z;
      p.z = depth;
      points_shrinked.push_back(p);
    }
}

#endif /* PCL_REGISTRATION_IMPL_ICP_LINES_HPP_ */
