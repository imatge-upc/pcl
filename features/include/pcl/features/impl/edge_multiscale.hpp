/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2012, Willow Garage, Inc.
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
 * $Id$
 */

#ifndef PCL_FEATURES_EDGE_MULTISCALE_H_
#define PCL_FEATURES_EDGE_MULTISCALE_H_

#include <utility>
#include <pcl/features/edge_multiscale.h>

//////////////////////////////////////////////////////////////////////////////////////////////
// Compute a local Reference Frame for a 3D feature; the output is stored in the "rf" matrix
template<typename PointInT, typename PointOutT, typename PointInternalT> void
pcl::MultiscaleEdgeDetector<PointInT, PointOutT, PointInternalT>::estimate (pcl::PointCloud<PointOutT> &output)
{
  if (!_cloud_set) {}
  if (!_kdtree_set) {computeKdtree();}
  typename pcl::PointCloud<PointInternalT>::Ptr neighbours (new pcl::PointCloud<PointInternalT>());
  typename pcl::ExtractIndices<PointInternalT> extract;
  extract.setInputCloud (_cloud);
  extract.setNegative (false);

  Eigen::RowVector3f mean, mean_step;
  Eigen::RowVector3f mean_intens, mean_intens_inv;
  Eigen::Vector3f eigenvals;
  Eigen::Matrix3f cov, cov_step;
  Eigen::Matrix3f M0, M1, M2, M3;
  float c3, c2, c1, c0;
  float sum_eigenvals, eigen0, eigen1;
  float pot3, pot2;
  float num, den;
  float prev_eigen;
  M0 = Eigen::Matrix3f::Zero();
  c3 = 1;
  M1 = Eigen::Matrix3f::Identity();
  Eigen::Matrix<float, Eigen::Dynamic, 3> centered;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig;

  boost::shared_ptr<std::vector<int> > pointIdxNKNSearch(new std::vector<int>(_max_K));
  boost::shared_ptr<std::vector<float> > pointNKNSquaredDistance(new std::vector<float>(_max_K));

  int iters = 0, iters_color = 0;
  float weight_step = 0, weight_color = 0;
  int previous_K=_min_K;
  Eigen::Matrix<float, Eigen::Dynamic, 3> neigh_matrix;
  Eigen::Matrix<float, Eigen::Dynamic, 3> intensities, intensities_inv;

  neigh_matrix = Eigen::Matrix<float, Eigen::Dynamic, 3>::Zero(_max_K,3);
  intensities = Eigen::Matrix<float, Eigen::Dynamic, 3>::Zero(_max_K,3);
  intensities_inv = Eigen::Matrix<float, Eigen::Dynamic, 3>::Zero(_max_K,3);

  int h, K, k;
  float stop_condition;
  int iters_nostep, iters_nocolor;  float w_step, w_color;
  float mean_distance, mean_distance_step;
  bool sum=false;
  float intens_sum, intens_sum_inv;
  std::cout << "Thresholds: " << _weight_level_thresh << " " << _weight_level_color_thresh << std::endl;
  std::cout << "K's: " << _min_K << " " << _max_K <<  " " << _max_K_color <<  std::endl;

  for (int j = 0; j < _cloud->points.size(); j++) {
    _kdtree.nearestKSearch(_cloud->points[j], _max_K, *pointIdxNKNSearch, *pointNKNSquaredDistance);
    Eigen::Map<Eigen::VectorXf> squaredistances(pointNKNSquaredDistance->data(), _max_K);
    extract.setIndices (pointIdxNKNSearch);
    extract.filter(*neighbours);
    iters = 0;
    iters_color = 0;
    weight_step = 0;
    weight_color = 0;
    previous_K = _min_K;
    iters_nocolor = 0;
    iters_nostep = 0;
    intens_sum = 0;
    intens_sum_inv = 0;

    // First iteration
    for (int h=0; h<_min_K; h++){ // Initialize first K values
        neigh_matrix(h,0) = neighbours->points[h].x;
        neigh_matrix(h,1) = neighbours->points[h].y;
        neigh_matrix(h,2) = neighbours->points[h].z;
        intensities(h,0) = neighbours->points[h].intensity;
        intensities(h,1) = neighbours->points[h].intensity;
        intensities(h,2) = neighbours->points[h].intensity;
        intensities_inv(h,0) = 1-neighbours->points[h].intensity;
        intensities_inv(h,1) = 1-neighbours->points[h].intensity;
        intensities_inv(h,2) = 1-neighbours->points[h].intensity;
        intens_sum += neighbours->points[h].intensity;
        intens_sum_inv += 1-neighbours->points[h].intensity;
        // std::cout << intens_sum << std::endl;
    }

    
    mean = neigh_matrix.block(0,0, _min_K, 3).colwise().mean();
    mean_distance = (squaredistances.block(0,0, _min_K, 1).colwise().mean())(0);
    mean_intens = neigh_matrix.block(0,0, _min_K, 3).cwiseProduct(
                    intensities.block(0,0, _min_K, 3)).colwise().sum()/intens_sum;
    mean_intens_inv = neigh_matrix.block(0,0, _min_K, 3).cwiseProduct(
                    intensities_inv.block(0,0, _min_K, 3)).colwise().sum()/intens_sum_inv;
    // w_color = (mean-mean_intens).squaredNorm()/(mean_distance/2);
    w_color = std::min((mean-mean_intens).squaredNorm(),(mean-mean_intens_inv).squaredNorm())/(mean_distance*0.36);
    centered = neigh_matrix.block(0,0, _min_K, 3).rowwise() - mean;
    // std::cout << "centered: " << std::endl << centered << std::endl;
    cov = (centered.transpose() * centered)/ float(_min_K);
    eig.compute(cov, Eigen::EigenvaluesOnly);
    eigenvals = eig.eigenvalues();
    sum_eigenvals = eigenvals.sum();
    eigen0 = eigenvals[0];
    w_step = eigen0/eigenvals.sum();

    if (w_step > _weight_level_thresh) {
        weight_step += 1.0;
        iters_nostep = 0;
    } else {
        iters_nostep++;
    }  
    if (w_color > _weight_level_color_thresh) {
        weight_color += 1.0;
        iters_nocolor = 0;
    } else {
        iters_nocolor++;
    }
    iters ++;
    iters_color ++;
    if (_single_k){
        PointOutT p2;
        p2.x = _cloud->points[j].x; p2.y = _cloud->points[j].y; p2.z = _cloud->points[j].z;
        p2.intensity = w_color*10;
        output.push_back(p2);

    } else {
        for (K=_min_K+_step_K; K<_max_K; K+=_step_K){
            if ( iters_nostep != iters || iters_nostep < _iters_stop_nostep || 
                 iters_nocolor != iters || iters_nocolor < _iters_stop_nocolor){
                for (h=previous_K; h<K; h++){
                    neigh_matrix(h,0) = neighbours->points[h].x;
                    neigh_matrix(h,1) = neighbours->points[h].y;
                    neigh_matrix(h,2) = neighbours->points[h].z;
                    intensities(h,0) = neighbours->points[h].intensity;
                    intensities(h,1) = neighbours->points[h].intensity;
                    intensities(h,2) = neighbours->points[h].intensity;
                    intensities_inv(h,0) = 1-neighbours->points[h].intensity;
                    intensities_inv(h,1) = 1-neighbours->points[h].intensity;
                    intensities_inv(h,2) = 1-neighbours->points[h].intensity;                        
                    intens_sum += neighbours->points[h].intensity;
                    intens_sum_inv += 1-neighbours->points[h].intensity;
                    // std::cout << intens_sum << std::endl;
                }
                // std::cout << "updated matrices" << std::endl;
                mean_step = neigh_matrix.block(previous_K,0, _step_K, 3).colwise().mean();
                mean = (previous_K*mean+_step_K*mean_step)/K;
                // mean = neigh_matrix.block(0,0, K, 3).colwise().mean();
                if (( iters_nocolor != iters || iters_nocolor < _iters_stop_nocolor) 
                    && K <= _max_K_color){
                    mean_distance_step = (squaredistances.block(previous_K,0, _step_K, 1).colwise().mean())(0);
                    mean_distance = (previous_K*mean_distance+_step_K*mean_distance_step)/K;
                    mean_intens = neigh_matrix.block(0,0, K, 3).cwiseProduct(
                                    intensities.block(0,0, K, 3)).colwise().sum()/intens_sum;
                    mean_intens_inv = neigh_matrix.block(0,0, K, 3).cwiseProduct(
                                    intensities_inv.block(0,0, K, 3)).colwise().sum()/intens_sum_inv;
                   w_color = std::min((mean-mean_intens).squaredNorm(),(mean-mean_intens_inv).squaredNorm())/(mean_distance*0.36);
                } else {w_color=0;}
                if (iters_nostep != iters || iters_nostep < _iters_stop_nostep){
                    centered = neigh_matrix.block(previous_K,0, _step_K, 3).rowwise() - mean;
                    cov_step = centered.transpose() * centered; // *Step_k
                    cov = (previous_K*cov+cov_step)/K;
                    // ToDO: update this!
                    centered = neigh_matrix.block(0,0, K, 3).rowwise() - mean;
                    cov = (centered.transpose() * centered) / float(K); // *Step_k

                    if (K>5){ // Suposarem que si ja tinc 5 punts no tinc una línia. Si la tinc, cagada.
                        c2 = -cov.trace();
                        M2 = cov*M1+c2*Eigen::Matrix3f::Identity();
                        c1 = -0.5*(cov*M2).trace();
                        M3 = cov*M2+c1*Eigen::Matrix3f::Identity();
                        c0 = -(1.0/3.0)*(cov*M3).trace();
                        sum_eigenvals = -c2;
                        stop_condition = abs(sum_eigenvals/1000000);
                        bool failed = false;
                        do {
                            prev_eigen = eigen0;
                            pot3 = pow(eigen0, 3.0); pot2 = pow(eigen0, 2.0);
                            num = c3*pot3+c2*pot2+c1*eigen0+c0;
                            den = 3*c3*pot2+2*c2*eigen0+c1;
                            eigen0 -= num/den;
                        } while(abs(prev_eigen-eigen0)> stop_condition);
                    } else {
                        eig.compute(cov, Eigen::EigenvaluesOnly);
                        eigenvals = eig.eigenvalues();
                        sum_eigenvals = eigenvals.sum();
                        eigen0 = eigenvals[0];
                        iters_nostep = 0;
                    }
                    w_step = eigen0/sum_eigenvals;
                } else {w_step=0;}
            } else {w_color=0; w_step=0;}
            sum=false;

            if (w_step > _weight_level_thresh) {
                weight_step += 1.0;
                iters_nostep = 0;
            } else {
                iters_nostep++;
            }  
            if (w_color > _weight_level_color_thresh) {
                weight_color += 1.0;
                iters_nocolor = 0;
            } else {
                iters_nocolor++;
            }
            iters++;
            if (K <= _max_K_color) iters_color ++;
            previous_K = K;
        }
       
        weight_step /= iters;
        weight_color /= iters_color;
        PointOutT p2;
        p2.x = _cloud->points[j].x; p2.y = _cloud->points[j].y; p2.z = _cloud->points[j].z;
        p2.intensity = std::max(weight_step, weight_color)*100;
        // p2.intensity=_cloud->points[j].intensity;
        output.push_back(p2);
    }
  }
}


template<typename PointInT, typename PointOutT, typename PointInternalT> float
pcl::MultiscaleEdgeDetector<PointInT, PointOutT, PointInternalT>::estimatePoint (PointInternalT &point)
{
  if (!_cloud_set) {}
  if (!_kdtree_set) {computeKdtree();}
  typename pcl::PointCloud<PointInternalT>::Ptr neighbours (new pcl::PointCloud<PointInternalT>());
  typename pcl::ExtractIndices<PointInternalT> extract;
  extract.setInputCloud (_cloud);
  extract.setNegative (false);

  Eigen::RowVector3f mean, mean_step;
  Eigen::RowVector3f mean_intens, mean_intens_inv;
  Eigen::Vector3f eigenvals;
  Eigen::Matrix3f cov, cov_step;
  Eigen::Matrix3f M0, M1, M2, M3;
  float c3, c2, c1, c0;
  float sum_eigenvals, eigen0, eigen1;
  float pot3, pot2;
  float num, den;
  float prev_eigen;
  M0 = Eigen::Matrix3f::Zero();
  c3 = 1;
  M1 = Eigen::Matrix3f::Identity();
  Eigen::Matrix<float, Eigen::Dynamic, 3> centered;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig;

  boost::shared_ptr<std::vector<int> > pointIdxNKNSearch(new std::vector<int>(_max_K));
  boost::shared_ptr<std::vector<float> > pointNKNSquaredDistance(new std::vector<float>(_max_K));

  int iters = 0, iters_color = 0;
  float weight_step = 0, weight_color = 0;
  int previous_K=_min_K;
  Eigen::Matrix<float, Eigen::Dynamic, 3> neigh_matrix;
  Eigen::Matrix<float, Eigen::Dynamic, 3> intensities, intensities_inv;

  neigh_matrix = Eigen::Matrix<float, Eigen::Dynamic, 3>::Zero(_max_K,3);
  intensities = Eigen::Matrix<float, Eigen::Dynamic, 3>::Zero(_max_K,3);
  intensities_inv = Eigen::Matrix<float, Eigen::Dynamic, 3>::Zero(_max_K,3);

  int h, K, k;
  float stop_condition;
  int iters_nostep, iters_nocolor;  float w_step, w_color;
  float mean_distance, mean_distance_step;
  bool sum=false;
  float intens_sum, intens_sum_inv;
    _kdtree.nearestKSearch(point, _max_K, *pointIdxNKNSearch, *pointNKNSquaredDistance);
    Eigen::Map<Eigen::VectorXf> squaredistances(pointNKNSquaredDistance->data(), _max_K);
    extract.setIndices (pointIdxNKNSearch);
    extract.filter(*neighbours);
    iters = 0;
    iters_color = 0;
    weight_step = 0;
    weight_color = 0;
    previous_K = _min_K;
    iters_nocolor = 0;
    iters_nostep = 0;
    intens_sum = 0;
    intens_sum_inv = 0;

    // First iteration
    for (int h=0; h<_min_K; h++){ // Initialize first K values
        neigh_matrix(h,0) = neighbours->points[h].x;
        neigh_matrix(h,1) = neighbours->points[h].y;
        neigh_matrix(h,2) = neighbours->points[h].z;
        intensities(h,0) = neighbours->points[h].intensity;
        intensities(h,1) = neighbours->points[h].intensity;
        intensities(h,2) = neighbours->points[h].intensity;
        intensities_inv(h,0) = 1-neighbours->points[h].intensity;
        intensities_inv(h,1) = 1-neighbours->points[h].intensity;
        intensities_inv(h,2) = 1-neighbours->points[h].intensity;
        intens_sum += neighbours->points[h].intensity;
        intens_sum_inv += 1-neighbours->points[h].intensity;
    }

    
    mean = neigh_matrix.block(0,0, _min_K, 3).colwise().mean();
    mean_distance = (squaredistances.block(0,0, _min_K, 1).colwise().mean())(0);
    mean_intens = neigh_matrix.block(0,0, _min_K, 3).cwiseProduct(
                    intensities.block(0,0, _min_K, 3)).colwise().sum()/intens_sum;
    mean_intens_inv = neigh_matrix.block(0,0, _min_K, 3).cwiseProduct(
                    intensities_inv.block(0,0, _min_K, 3)).colwise().sum()/intens_sum_inv;
    // w_color = (mean-mean_intens).squaredNorm()/(mean_distance/2);
    w_color = std::min((mean-mean_intens).squaredNorm(),(mean-mean_intens_inv).squaredNorm())/(mean_distance*0.36);
    centered = neigh_matrix.block(0,0, _min_K, 3).rowwise() - mean;
    // std::cout << "centered: " << std::endl << centered << std::endl;
    cov = (centered.transpose() * centered)/ float(_min_K);
    eig.compute(cov, Eigen::EigenvaluesOnly);
    eigenvals = eig.eigenvalues();
    sum_eigenvals = eigenvals.sum();
    eigen0 = eigenvals[0];
    w_step = eigen0/eigenvals.sum();

    if (w_step > _weight_level_thresh) {
        weight_step += 1.0;
        iters_nostep = 0;
    } else {
        iters_nostep++;
    }  
    if (w_color > _weight_level_color_thresh) {
        weight_color += 1.0;
        iters_nocolor = 0;
    } else {
        iters_nocolor++;
    }
    iters ++;
    iters_color ++;
    if (_single_k){
      return std::max(w_color, w_step);

    } else {
      for (K=_min_K+_step_K; K<_max_K; K+=_step_K){
          if ( iters_nostep != iters || iters_nostep < _iters_stop_nostep || 
               iters_nocolor != iters || iters_nocolor < _iters_stop_nocolor){
              for (h=previous_K; h<K; h++){
                  neigh_matrix(h,0) = neighbours->points[h].x;
                  neigh_matrix(h,1) = neighbours->points[h].y;
                  neigh_matrix(h,2) = neighbours->points[h].z;
                  intensities(h,0) = neighbours->points[h].intensity;
                  intensities(h,1) = neighbours->points[h].intensity;
                  intensities(h,2) = neighbours->points[h].intensity;
                  intensities_inv(h,0) = 1-neighbours->points[h].intensity;
                  intensities_inv(h,1) = 1-neighbours->points[h].intensity;
                  intensities_inv(h,2) = 1-neighbours->points[h].intensity;                        
                  intens_sum += neighbours->points[h].intensity;
                  intens_sum_inv += 1-neighbours->points[h].intensity;
              }
              // std::cout << "updated matrices" << std::endl;
              mean_step = neigh_matrix.block(previous_K,0, _step_K, 3).colwise().mean();
              mean = (previous_K*mean+_step_K*mean_step)/K;
              // mean = neigh_matrix.block(0,0, K, 3).colwise().mean();
              if (( iters_nocolor != iters || iters_nocolor < _iters_stop_nocolor) 
                  && K <= _max_K_color){
                  mean_distance_step = (squaredistances.block(previous_K,0, _step_K, 1).colwise().mean())(0);
                  mean_distance = (previous_K*mean_distance+_step_K*mean_distance_step)/K;
                  mean_intens = neigh_matrix.block(0,0, K, 3).cwiseProduct(
                                  intensities.block(0,0, K, 3)).colwise().sum()/intens_sum;
                  mean_intens_inv = neigh_matrix.block(0,0, K, 3).cwiseProduct(
                                  intensities_inv.block(0,0, K, 3)).colwise().sum()/intens_sum_inv;
                 w_color = std::min((mean-mean_intens).squaredNorm(),(mean-mean_intens_inv).squaredNorm())/(mean_distance*0.36);
              } else {w_color=0;}
              if (iters_nostep != iters || iters_nostep < _iters_stop_nostep){
                  centered = neigh_matrix.block(previous_K,0, _step_K, 3).rowwise() - mean;
                  cov_step = centered.transpose() * centered; // *Step_k
                  cov = (previous_K*cov+cov_step)/K;
                  // ToDO: update this!
                  centered = neigh_matrix.block(0,0, K, 3).rowwise() - mean;
                  cov = (centered.transpose() * centered) / float(K); // *Step_k

                  if (K>5){ // Suposarem que si ja tinc 5 punts no tinc una línia. Si la tinc, cagada.
                      c2 = -cov.trace();
                      M2 = cov*M1+c2*Eigen::Matrix3f::Identity();
                      c1 = -0.5*(cov*M2).trace();
                      M3 = cov*M2+c1*Eigen::Matrix3f::Identity();
                      c0 = -(1.0/3.0)*(cov*M3).trace();
                      sum_eigenvals = -c2;
                      stop_condition = abs(sum_eigenvals/1000000);
                      bool failed = false;
                      do {
                          prev_eigen = eigen0;
                          pot3 = pow(eigen0, 3.0); pot2 = pow(eigen0, 2.0);
                          num = c3*pot3+c2*pot2+c1*eigen0+c0;
                          den = 3*c3*pot2+2*c2*eigen0+c1;
                          eigen0 -= num/den;
                      } while(abs(prev_eigen-eigen0)> stop_condition);
                  } else {
                      eig.compute(cov, Eigen::EigenvaluesOnly);
                      eigenvals = eig.eigenvalues();
                      sum_eigenvals = eigenvals.sum();
                      eigen0 = eigenvals[0];
                      iters_nostep = 0;
                  }
                  w_step = eigen0/sum_eigenvals;
              } else {w_step=0;}
          } else {w_color=0; w_step=0;}
          sum=false;

          if (w_step > _weight_level_thresh) {
              weight_step += 1.0;
              iters_nostep = 0;
          } else {
              iters_nostep++;
          }  
          if (w_color > _weight_level_color_thresh) {
              weight_color += 1.0;
              iters_nocolor = 0;
          } else {
              iters_nocolor++;
          }
          iters++;
          if (K <= _max_K_color) iters_color ++;
          previous_K = K;
      }
      
      weight_step /= iters;
      weight_color /= iters_color;
      return std::max(weight_step, weight_color);
  }
}

template<typename PointInT, typename PointOutT, typename PointInternalT> void
pcl::MultiscaleEdgeDetector<PointInT, PointOutT, PointInternalT>::setInputCloud (pcl::PointCloud<PointInT> &cloud)
{
  _cloud = typename pcl::PointCloud<PointInternalT>::Ptr(new pcl::PointCloud<PointInternalT>);
  for (int i=0; i< cloud.points.size(); i++){
      PointInternalT p;
      p.x = cloud.points[i].x;
      p.y = cloud.points[i].y;
      p.z = cloud.points[i].z;
      p.intensity = (0.299*float(cloud.points[i].r) +
                     0.587*float(cloud.points[i].g) +
                     0.114*float(cloud.points[i].b))/255.0;
      _cloud->push_back(p);
  }
  _cloud_set = true;
}

// template<typename PointInT, typename PointOutT, typename PointInternalT> void
// pcl::MultiscaleEdgeDetector<PointInT, PointOutT, PointInternalT>::setInputCloud (pcl::PointCloud<PointInT> &cloud)
// {
//   _cloud = pcl::PointCloud<PointInternalT>::Ptr(cloud);
// }



#define PCL_INSTANTIATE_MultiscaleEdgeDetector(T,IntT,OutT) template class PCL_EXPORTS pcl::MultiscaleEdgeDetector<T,IntT, OutT>;

#endif    // PCL_FEATURES_EDGE_MULTISCALE_H_

