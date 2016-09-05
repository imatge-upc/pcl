/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
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
 */


#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/io/pcd_io.h>
 #include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/common/transforms.h>
#include <vtkVersion.h>
#include <vtkPLYReader.h>
#include <vtkOBJReader.h>
//#include <vtkOBJImporter.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/console/print.h>
#include <pcl/console/parse.h>
#include <pcl/kdtree/kdtree_flann.h>

 #include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/uniform.hpp>

inline double
uniform_deviate (int seed)
{
  double ran = seed * (1.0 / (RAND_MAX + 1.0));
  return ran;
}

inline void
randomPointTriangle (float a1, float a2, float a3, float b1, float b2, float b3, float c1, float c2, float c3,
                     Eigen::Vector4f& p)
{
  float r1 = static_cast<float> (uniform_deviate (rand ()));
  float r2 = static_cast<float> (uniform_deviate (rand ()));
  float r1sqr = sqrtf (r1);
  float OneMinR1Sqr = (1 - r1sqr);
  float OneMinR2 = (1 - r2);
  a1 *= OneMinR1Sqr;
  a2 *= OneMinR1Sqr;
  a3 *= OneMinR1Sqr;
  b1 *= OneMinR2;
  b2 *= OneMinR2;
  b3 *= OneMinR2;
  c1 = r1sqr * (r2 * c1 + b1) + a1;
  c2 = r1sqr * (r2 * c2 + b2) + a2;
  c3 = r1sqr * (r2 * c3 + b3) + a3;
  p[0] = c1;
  p[1] = c2;
  p[2] = c3;
  p[3] = 0;
}

inline void
randPSurface (vtkPolyData * polydata, std::vector<double> * cumulativeAreas, double totalArea, Eigen::Vector4f& p, Eigen::Vector3f& c,
    pcl::PointCloud<pcl::PointXYZRGB> & color_cloud, pcl::KdTreeFLANN<pcl::PointXYZRGB> &kdtree)
{
  float r = static_cast<float> (uniform_deviate (rand ()) * totalArea);

  std::vector<double>::iterator low = std::lower_bound (cumulativeAreas->begin (), cumulativeAreas->end (), r);
  vtkIdType el = vtkIdType (low - cumulativeAreas->begin ());

  double A[3], B[3], C[3];
  vtkIdType npts = 0;
  vtkIdType *ptIds = NULL;
  polydata->GetCellPoints (el, npts, ptIds);
  polydata->GetPoint (ptIds[0], A);
  polydata->GetPoint (ptIds[1], B);
  polydata->GetPoint (ptIds[2], C);
  randomPointTriangle (float (A[0]), float (A[1]), float (A[2]), 
                       float (B[0]), float (B[1]), float (B[2]), 
                       float (C[0]), float (C[1]), float (C[2]), p);
  float d[3], sum_d;
  d[0] = sqrt((float(A[0])-c[0])*(float(A[0])-c[0])+(float(A[1])-c[1])*(float(A[1])-c[1])+(float(A[2])-c[2])*(float(A[2])-c[2]));
  d[1] = sqrt((float(B[0])-c[0])*(float(B[0])-c[0])+(float(B[1])-c[1])*(float(B[1])-c[1])+(float(B[2])-c[2])*(float(B[2])-c[2]));
  d[2] = sqrt((float(C[0])-c[0])*(float(C[0])-c[0])+(float(C[1])-c[1])*(float(C[1])-c[1])+(float(C[2])-c[2])*(float(C[2])-c[2]));
  sum_d = d[0]+d[1]+d[2];
  d[0] /= sum_d; d[1] /=sum_d; d[2] /=sum_d;
  c[0] = 0; c[1] = 0; c[2] = 0;
  c[0] += color_cloud.points[ptIds[0]].r*d[0];
  c[1] += color_cloud.points[ptIds[0]].g*d[0];
  c[2] += color_cloud.points[ptIds[0]].b*d[0];  
  c[0] += color_cloud.points[ptIds[1]].r*d[1];
  c[1] += color_cloud.points[ptIds[1]].g*d[1];
  c[2] += color_cloud.points[ptIds[1]].b*d[1];  
  c[0] += color_cloud.points[ptIds[2]].r*d[2];
  c[1] += color_cloud.points[ptIds[2]].g*d[2];
  c[2] += color_cloud.points[ptIds[2]].b*d[2]; 

}

void
uniform_sampling (vtkSmartPointer<vtkPolyData> polydata, size_t n_samples, 
  pcl::PointCloud<pcl::PointXYZRGB> & color_cloud,
  pcl::PointCloud<pcl::PointXYZRGB> & cloud_out,
  float color_noise, float normal_noise)
{
  pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
  kdtree.setInputCloud(typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr(&color_cloud));
  polydata->BuildCells ();
  vtkSmartPointer<vtkCellArray> cells = polydata->GetPolys ();

  double p1[3], p2[3], p3[3], totalArea = 0;
  std::vector<double> cumulativeAreas (cells->GetNumberOfCells (), 0);
  size_t i = 0;
  vtkIdType npts = 0, *ptIds = NULL;
  for (cells->InitTraversal (); cells->GetNextCell (npts, ptIds); i++)
  {
    polydata->GetPoint (ptIds[0], p1);
    polydata->GetPoint (ptIds[1], p2);
    polydata->GetPoint (ptIds[2], p3);
    totalArea += vtkTriangle::TriangleArea (p1, p2, p3);
    cumulativeAreas[i] = totalArea;
  }

  cloud_out.points.resize (n_samples);
  cloud_out.width = static_cast<pcl::uint32_t> (n_samples);
  cloud_out.height = 1;


  boost::mt19937 *rng = new boost::mt19937();
  rng->seed(time(NULL));

  boost::normal_distribution<> distribution_n(0, normal_noise);
  boost::variate_generator< boost::mt19937, boost::normal_distribution<> > dist_n(*rng, distribution_n);
    
  for (i = 0; i < n_samples; i++)
  {
    Eigen::Vector4f p;
    Eigen::Vector3f c;
    randPSurface (polydata, &cumulativeAreas, totalArea, p, c,color_cloud, kdtree);
    float d_c = 0;
    if (color_noise>0) d_c = (float(rand())*2*color_noise/float(RAND_MAX)-color_noise);
    float c0 = c[0], c1 = c[1], c2=c[2];
    // std::cout << c0 << std::endl;
    // std::cout <<d_c << std::endl;

    // std::cout << int(std::max(0.0f, std::min(255.0f, c0 + d_c))+0.5f) << std::endl << std::endl;
    cloud_out.points[i].x = p[0] + dist_n();
    cloud_out.points[i].y = p[1] + dist_n();
    cloud_out.points[i].z = p[2] + dist_n();
    cloud_out.points[i].r = int(std::max(0.0f, std::min(255.0f, c0 + d_c))+0.5f);
    cloud_out.points[i].g = int(std::max(0.0f, std::min(255.0f, c1 + d_c))+0.5f);
    cloud_out.points[i].b = int(std::max(0.0f, std::min(255.0f, c2 + d_c))+0.5f);
  }
}

using namespace pcl;
using namespace pcl::io;
using namespace pcl::console;

int default_number_samples = 100000;
float default_leaf_size = 0.01f;

void
printHelp (int, char **argv)
{
  print_error ("Syntax is: %s input.{ply,obj} output.pcd <options>\n", argv[0]);
  print_info ("  where options are:\n");
  print_info ("                     -n_samples X      = number of samples (default: ");
  print_value ("%d", default_number_samples);
  print_info (")\n");
  print_info (
              "                     -leaf_size X  = the XYZ leaf size for the VoxelGrid -- for data reduction (default: ");
  print_value ("%f", default_leaf_size);
  print_info (" m)\n");
}

/* ---[ */
int
main (int argc, char **argv)
{
  print_info ("Convert a CAD model to a point cloud using uniform sampling. For more information, use: %s -h\n",
              argv[0]);

  if (argc < 3)
  {
    printHelp (argc, argv);
    return (-1);
  }

  // Parse command line arguments
  int SAMPLE_POINTS_ = default_number_samples;
  parse_argument (argc, argv, "--n_samples", SAMPLE_POINTS_);
  float leaf_size = default_leaf_size;
  parse_argument (argc, argv, "--leaf_size", leaf_size);  
  float color_noise = 0;
  float normal_noise =0;
  parse_argument (argc, argv, "--color_noise", color_noise);
  parse_argument (argc, argv, "--normal_noise", normal_noise);

  // Parse the command line arguments for .ply and PCD files
  std::vector<int> pcd_file_indices = parse_file_extension_argument (argc, argv, ".pcd");
  if (pcd_file_indices.size () != 1)
  {
    print_error ("Need a single output PCD file to continue.\n");
    return (-1);
  }
  std::vector<int> ply_file_indices = parse_file_extension_argument (argc, argv, ".ply");
  std::vector<int> obj_file_indices = parse_file_extension_argument (argc, argv, ".obj");
  if (ply_file_indices.size () != 1 && obj_file_indices.size () != 1)
  {
    print_error ("Need a single input PLY/OBJ file to continue.\n");
    return (-1);
  }
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr points_color(new pcl::PointCloud<pcl::PointXYZRGB>);
  vtkSmartPointer<vtkPolyData> polydata1 = vtkSmartPointer<vtkPolyData>::New ();;
  if (ply_file_indices.size () == 1)
  {
    pcl::PolygonMesh mesh;
    pcl::io::loadPolygonFilePLY (argv[ply_file_indices[0]], mesh);
    pcl::io::mesh2vtk (mesh, polydata1);
    // std::cout << mesh << std::endl;
    pcl::io::loadPLYFile<pcl::PointXYZRGB> (argv[ply_file_indices[0]], *points_color);
  }
  else if (obj_file_indices.size () == 1)
  {
    vtkSmartPointer<vtkOBJReader> readerQuery = vtkSmartPointer<vtkOBJReader>::New ();
    readerQuery->SetFileName (argv[obj_file_indices[0]]);
    readerQuery->Update ();
    polydata1 = readerQuery->GetOutput ();
  }

  //make sure that the polygons are triangles!
  vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New ();
#if VTK_MAJOR_VERSION < 6
  triangleFilter->SetInput (polydata1);
#else
  triangleFilter->SetInputData (polydata1);
#endif
  triangleFilter->Update ();

  vtkSmartPointer<vtkPolyDataMapper> triangleMapper = vtkSmartPointer<vtkPolyDataMapper>::New ();
  triangleMapper->SetInputConnection (triangleFilter->GetOutputPort ());
  triangleMapper->Update();
  polydata1 = triangleMapper->GetInput();

  bool INTER_VIS = false;
  bool VIS = false;

  if (INTER_VIS)
  {
    visualization::PCLVisualizer vis;
    vis.addModelFromPolyData (polydata1, "mesh1", 0);
    vis.setRepresentationToSurfaceForAllActors ();
    vis.spin();
  }

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_1 (new pcl::PointCloud<pcl::PointXYZRGB>);
  uniform_sampling (polydata1, SAMPLE_POINTS_, *points_color, *cloud_1, color_noise, normal_noise);

  if (INTER_VIS)
  {
    visualization::PCLVisualizer vis_sampled;
    vis_sampled.addPointCloud (cloud_1);
    vis_sampled.spin ();
  }

  // Voxelgrid
  VoxelGrid<PointXYZRGB> grid_;
  grid_.setInputCloud (cloud_1);
  grid_.setLeafSize (leaf_size, leaf_size, leaf_size);

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr res(new pcl::PointCloud<pcl::PointXYZRGB>);
  grid_.filter (*res);

  if (VIS)
  {
    visualization::PCLVisualizer vis3 ("VOXELIZED SAMPLES CLOUD");
    vis3.addPointCloud (res);
    vis3.spin ();
  }

  savePCDFileASCII (argv[pcd_file_indices[0]], *res);
}
