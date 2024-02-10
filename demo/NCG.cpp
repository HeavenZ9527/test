#include <iostream> //标准输入输出流
#include <pcl/io/pcd_io.h> //PCL的PCD格式文件的输入输出头文件
#include <pcl/point_types.h> //PCL对各种格式的点的支持头文件
#include <pcl/visualization/cloud_viewer.h>//点云查看窗口头文件
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
//#include <pcl/filters/extract_indices.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/common/transforms.h>
#include <pcl/common/transformation_from_correspondences.h>
#include <pcl/console/time.h>
#include <pcl/common/transforms.h> 
#include <stdlib.h> 
#include <random>
#include <iomanip>//设置精度
#include <cstdlib>
#include <cmath>
#include <limits>
#include <time.h>

#define search_range 0.06 //聚类最大距离
#define min_cluster 30 //最小聚类大小
#define curv_thresh 5 //曲率阈值
#define reg_size 4 // 参与配准的点数
#define MIU 0     
#define SIGMA 0.5 
#define desc_type 3 //特征描述子类型
 

struct feature_desc
{
	int idx;
	float s1;
	float s2;
	float s3;
	float i1;
	float i2;
	float i3;
};
struct feature_desc_4
{
	int idx;
	float s1;
	float s2;
	float s3;
	float s4;
	float s5;
	float s6;
	float i1;
	float i2;
	float i3;
	float i4;
};

struct feature_desc_5
{
	int idx;
	float s1;
	float s2;
	float s3;
	float s4;
	float s5;
	float s6;
	float s7;
	float s8;
	float s9;
	float s10;
	float i1;
	float i2;
	float i3;
	float i4;
	float i5;
};

struct feature_desc_6
{
	int idx;
	float s1;
	float s2;
	float s3;
	float s4;
	float s5;
	float s6;
	float s7;
	float s8;
	float s9;
	float s10;
	float s11;
	float s12;
	float s13;
	float s14;
	float s15;
	float i1;
	float i2;
	float i3;
	float i4;
	float i5;
	float i6;
};

struct point_pair
{
	int idx_1;
	int idx_2;
	float d;
};
bool GreaterSort(point_pair a, point_pair b) { return (a.d>b.d); }
bool LessSort(point_pair a, point_pair b) { return (a.d<b.d); }

int filter_intensity(pcl::PointCloud<pcl::PointXYZI>::Ptr ptr_cloud, pcl::PointCloud<pcl::PointXYZI>::Ptr cloudOut, float thrsh)
{
	std::vector<int> indexs;
	int k = 0;
	for (size_t i = 0; i < ptr_cloud->points.size(); i++)
	{
		if (ptr_cloud->points[i].intensity>thrsh && ptr_cloud->points[i].intensity < thrsh*500)
		{
			indexs.push_back(i);
			k++;
		}

	}
	std::cout << " finished! with j= " << k << std::endl;

	pcl::copyPointCloud(*ptr_cloud, indexs, *cloudOut);

	return 0;
}
int cluster_mean(pcl::PointCloud<pcl::PointXYZI>::Ptr cloudOut, pcl::PointCloud<pcl::PointXYZI>::Ptr * add_cloud, pcl::PointCloud<pcl::PointXYZI>::Ptr * cloud_mean_points)
{
	pcl::PCDWriter writer;
	pcl::search::KdTree<pcl::PointXYZI>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZI>);
	tree->setInputCloud(cloudOut);

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZI> ec;   
	ec.setClusterTolerance(search_range);
	ec.setMinClusterSize(min_cluster);                
	ec.setMaxClusterSize(20000);              
	ec.setSearchMethod(tree);                   
	ec.setInputCloud(cloudOut);
	ec.extract(cluster_indices);           
										  
	float x_sum, y_sum, z_sum, intensity_sum = 0.0;
	int count = 0;

	for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
	{ 
		x_sum = 0.0;
		y_sum = 0.0;
		z_sum = 0.0;
		intensity_sum = 0.0;
		count = 0;
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_cluster(new pcl::PointCloud<pcl::PointXYZI>);
		for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); ++pit)
			//设置保存点云的属性问题
		{
			cloud_cluster->points.push_back(cloudOut->points[*pit]); //*
			x_sum += cloudOut->points[*pit].x;
			y_sum += cloudOut->points[*pit].y;
			z_sum += cloudOut->points[*pit].z;
			intensity_sum += cloudOut->points[*pit].intensity;
			count++;
		}

		pcl::PointXYZI mean_point;
		mean_point.x = x_sum / count;
		mean_point.y = y_sum / count;
		mean_point.z = z_sum / count;
		mean_point.intensity = intensity_sum / count;
		(*cloud_mean_points)->points.push_back(mean_point);

		cloud_cluster->width = cloud_cluster->points.size();
		cloud_cluster->height = 1;
		cloud_cluster->is_dense = true;
		
		std::cout << "PointCloud representing the Cluster: " << cloud_cluster->points.size() << " data points." << std::endl;
		**add_cloud += *cloud_cluster;
	}
	(*cloud_mean_points)->width = (*cloud_mean_points)->points.size();
	(*cloud_mean_points)->height = 1;
	(*cloud_mean_points)->is_dense = true;
	return 0;
}

int feature_discriptor(std::vector<feature_desc> &fea_vec, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mean_points)
{
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;               
	kdtree.setInputCloud(cloud_mean_points);                            
		int K = 3;
		std::vector<int> pointIdxNKNSearch(K);                    
		std::vector<float> pointNKNSquaredDistance(K);          

		for (size_t i = 0; i < cloud_mean_points->points.size(); ++i)
		{


			if (kdtree.nearestKSearch(cloud_mean_points->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)                                                             //执行k近邻搜索
			{
				feature_desc temp_fea;
				temp_fea.idx = i;
				temp_fea.s1 = sqrt(pointNKNSquaredDistance[1]);
				temp_fea.s2 = sqrt(pointNKNSquaredDistance[2]);
				temp_fea.s3 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[2]].x, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[2]].y, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[2]].z, 2));
				temp_fea.i1 = cloud_mean_points->points[i].intensity;
				temp_fea.i2 = cloud_mean_points->points[pointIdxNKNSearch[1]].intensity;
				temp_fea.i3 = cloud_mean_points->points[pointIdxNKNSearch[2]].intensity;

				fea_vec.push_back(temp_fea);

			}

		}
		//生成特征描述子
		return 0;

}

int feature_discriptor_4(std::vector<feature_desc_4> &fea_vec, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mean_points)
{
	//生成特征描述子
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;              
	kdtree.setInputCloud(cloud_mean_points);                          
		int K = 4;
		std::vector<int> pointIdxNKNSearch(K);                   
		std::vector<float> pointNKNSquaredDistance(K);           

		for (size_t i = 0; i < cloud_mean_points->points.size(); ++i)
		{

			if (kdtree.nearestKSearch(cloud_mean_points->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)                                                             //执行k近邻搜索
			{
				feature_desc_4 temp_fea;
				temp_fea.idx = i;
				temp_fea.s1 = sqrt(pointNKNSquaredDistance[1]);
				temp_fea.s2 = sqrt(pointNKNSquaredDistance[2]);
				temp_fea.s3 = sqrt(pointNKNSquaredDistance[3]);
				temp_fea.s4 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[2]].x, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[2]].y, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[2]].z, 2));
				temp_fea.s5 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[3]].x, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[3]].y, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[3]].z, 2));
				temp_fea.s6 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[2]].x - cloud_mean_points->points[pointIdxNKNSearch[3]].x, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[2]].y - cloud_mean_points->points[pointIdxNKNSearch[3]].y, 2) +
					pow(cloud_mean_points->points[pointIdxNKNSearch[2]].z - cloud_mean_points->points[pointIdxNKNSearch[3]].z, 2));

				fea_vec.push_back(temp_fea);

			}

		}

		//生成特征描述子
		return 0;
}

int feature_discriptor_5(std::vector<feature_desc_5> &fea_vec, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mean_points)
{
	
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;              
	kdtree.setInputCloud(cloud_mean_points);                            
																	
	int K = 5;
	std::vector<int> pointIdxNKNSearch(K);                  
	std::vector<float> pointNKNSquaredDistance(K);     
	for (size_t i = 0; i < cloud_mean_points->points.size(); ++i)
	{
		if (kdtree.nearestKSearch(cloud_mean_points->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)                                                             //执行k近邻搜索
		{
			feature_desc_5 temp_fea;
			temp_fea.idx = i;
			temp_fea.s1 = sqrt(pointNKNSquaredDistance[1]);
			temp_fea.s2 = sqrt(pointNKNSquaredDistance[2]);
			temp_fea.s3 = sqrt(pointNKNSquaredDistance[3]);
			temp_fea.s4 = sqrt(pointNKNSquaredDistance[4]);
			temp_fea.s5 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[2]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[2]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[2]].z, 2));
			temp_fea.s6 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[3]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[3]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[3]].z, 2));
			temp_fea.s7 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[4]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[4]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[4]].z, 2));
			temp_fea.s8 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[2]].x - cloud_mean_points->points[pointIdxNKNSearch[3]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].y - cloud_mean_points->points[pointIdxNKNSearch[3]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].z - cloud_mean_points->points[pointIdxNKNSearch[3]].z, 2));
			temp_fea.s9 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[2]].x - cloud_mean_points->points[pointIdxNKNSearch[4]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].y - cloud_mean_points->points[pointIdxNKNSearch[4]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].z - cloud_mean_points->points[pointIdxNKNSearch[4]].z, 2));
			temp_fea.s10 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[3]].x - cloud_mean_points->points[pointIdxNKNSearch[4]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[3]].y - cloud_mean_points->points[pointIdxNKNSearch[4]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[3]].z - cloud_mean_points->points[pointIdxNKNSearch[4]].z, 2));

			fea_vec.push_back(temp_fea);

		}
	}
	//生成特征描述子
	return 0;
}

int feature_discriptor_6(std::vector<feature_desc_6> &fea_vec, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mean_points)
{
	//生成特征描述子
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;               
	kdtree.setInputCloud(cloud_mean_points);                          
																	
	int K = 6;
	std::vector<int> pointIdxNKNSearch(K);                    
	std::vector<float> pointNKNSquaredDistance(K);           

	for (size_t i = 0; i < cloud_mean_points->points.size(); ++i)
	{


		if (kdtree.nearestKSearch(cloud_mean_points->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)                                                             //执行k近邻搜索
		{
			feature_desc_6 temp_fea;
			temp_fea.idx = i;
			temp_fea.s1 = sqrt(pointNKNSquaredDistance[1]);
			temp_fea.s2 = sqrt(pointNKNSquaredDistance[2]);
			temp_fea.s3 = sqrt(pointNKNSquaredDistance[3]);
			temp_fea.s4 = sqrt(pointNKNSquaredDistance[4]);
			temp_fea.s5 = sqrt(pointNKNSquaredDistance[5]);
			temp_fea.s6 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[2]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[2]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[2]].z, 2));
			temp_fea.s7 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[3]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[3]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[3]].z, 2));
			temp_fea.s8 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[4]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[4]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[4]].z, 2));
			temp_fea.s9 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[1]].x - cloud_mean_points->points[pointIdxNKNSearch[5]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].y - cloud_mean_points->points[pointIdxNKNSearch[5]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[1]].z - cloud_mean_points->points[pointIdxNKNSearch[5]].z, 2));
			temp_fea.s10 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[2]].x - cloud_mean_points->points[pointIdxNKNSearch[3]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].y - cloud_mean_points->points[pointIdxNKNSearch[3]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].z - cloud_mean_points->points[pointIdxNKNSearch[3]].z, 2));
			temp_fea.s11 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[2]].x - cloud_mean_points->points[pointIdxNKNSearch[4]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].y - cloud_mean_points->points[pointIdxNKNSearch[4]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].z - cloud_mean_points->points[pointIdxNKNSearch[4]].z, 2));
			temp_fea.s12 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[2]].x - cloud_mean_points->points[pointIdxNKNSearch[5]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].y - cloud_mean_points->points[pointIdxNKNSearch[5]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[2]].z - cloud_mean_points->points[pointIdxNKNSearch[5]].z, 2));
			temp_fea.s13 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[3]].x - cloud_mean_points->points[pointIdxNKNSearch[4]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[3]].y - cloud_mean_points->points[pointIdxNKNSearch[4]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[3]].z - cloud_mean_points->points[pointIdxNKNSearch[4]].z, 2));
			temp_fea.s14 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[3]].x - cloud_mean_points->points[pointIdxNKNSearch[5]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[3]].y - cloud_mean_points->points[pointIdxNKNSearch[5]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[3]].z - cloud_mean_points->points[pointIdxNKNSearch[5]].z, 2));
			temp_fea.s15 = sqrt(pow(cloud_mean_points->points[pointIdxNKNSearch[4]].x - cloud_mean_points->points[pointIdxNKNSearch[5]].x, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[4]].y - cloud_mean_points->points[pointIdxNKNSearch[5]].y, 2) +
				pow(cloud_mean_points->points[pointIdxNKNSearch[4]].z - cloud_mean_points->points[pointIdxNKNSearch[5]].z, 2));

			fea_vec.push_back(temp_fea);
		}
	}
	//生成特征描述子
	return 0;
}



float tonimoto(feature_desc x1, feature_desc x2)
{
	cout << "x1 : " << x1.s1 << " " << x1.s2 << " " << x1.s3 << " " << x1.i1 << " " << x1.i2 << " " << x1.i3 << " " << endl;
	cout << "x2 : " << x2.s1 << " " << x2.s2 << " " << x2.s3 << " " << x2.i1 << " " << x2.i2 << " " << x2.i3 << " " << endl;
	float ab = x1.s1*x2.s1 + x1.s2*x2.s2 + x1.s3*x2.s3 + x1.i1*x2.i1 + x1.i2*x2.i2 + x1.i3*x2.i3;
	float a2 = x1.s1*x1.s1 + x1.s2*x1.s2 + x1.s3*x1.s3 + x1.i1*x1.i1 + x1.i2*x1.i2 + x1.i3*x1.i3;
	float b2 = x2.s1*x2.s1 + x2.s2*x2.s2 + x2.s3*x2.s3 + x2.i1*x2.i1 + x2.i2*x2.i2 + x2.i3*x2.i3;
	return  ab / (a2 + b2 - ab);

}

float euclidian(feature_desc x1, feature_desc x2)
{
	
	float ab = sqrt(pow((x1.s1 - x2.s1)/(x1.s1 < x2.s1 ? x1.s1 : x2.s1), 2)*pow(1000, abs(x1.s1 - x2.s1) / (x1.s1 < x2.s1 ? x1.s1 : x2.s1)) +
		pow((x1.s2 - x2.s2)/(x1.s2 < x2.s2 ? x1.s2 : x2.s2), 2)*pow(1000, abs(x1.s2 - x2.s2) / (x1.s2 < x2.s2 ? x1.s2 : x2.s2)) +
		pow((x1.s3 - x2.s3)/(x1.s3 < x2.s3 ? x1.s3 : x2.s3), 2)*pow(1000, abs(x1.s3 - x2.s3) / (x1.s3 < x2.s3 ? x1.s3 : x2.s3)));
		
	return  ab;
}

float euclidian_4(feature_desc_4 x1, feature_desc_4 x2)
{
	
	float ab = abs(x1.s1 - x2.s1) / (x1.s1 < x2.s1 ? x1.s1 : x2.s1) +
		abs(x1.s2 - x2.s2) / (x1.s2 < x2.s2 ? x1.s2 : x2.s2) +
		abs(x1.s3 - x2.s3) / (x1.s3 < x2.s3 ? x1.s3 : x2.s3)+ 
		abs(x1.s4 - x2.s4) / (x1.s4 < x2.s4 ? x1.s4 : x2.s4) +
		abs(x1.s5 - x2.s5) / (x1.s5 < x2.s5 ? x1.s5 : x2.s5) +
		abs(x1.s6 - x2.s6) / (x1.s6 < x2.s6 ? x1.s6 : x2.s6);

	return  ab;
}

float euclidian_5(feature_desc_5 x1, feature_desc_5 x2)
{
	
	float ab = abs(x1.s1 - x2.s1) / (x1.s1 < x2.s1 ? x1.s1 : x2.s1) +
		abs(x1.s2 - x2.s2) / (x1.s2 < x2.s2 ? x1.s2 : x2.s2) +
		abs(x1.s3 - x2.s3) / (x1.s3 < x2.s3 ? x1.s3 : x2.s3) +
		abs(x1.s4 - x2.s4) / (x1.s4 < x2.s4 ? x1.s4 : x2.s4) +
		abs(x1.s5 - x2.s5) / (x1.s5 < x2.s5 ? x1.s5 : x2.s5) +
		abs(x1.s6 - x2.s6) / (x1.s6 < x2.s6 ? x1.s6 : x2.s6) +
		abs(x1.s7 - x2.s7) / (x1.s7 < x2.s7 ? x1.s7 : x2.s7) +
		abs(x1.s8 - x2.s8) / (x1.s8 < x2.s8 ? x1.s8 : x2.s8) +
		abs(x1.s9 - x2.s9) / (x1.s9 < x2.s9 ? x1.s9 : x2.s9) +
		abs(x1.s10 - x2.s10) / (x1.s10 < x2.s10 ? x1.s10 : x2.s10);

	return  ab;
}

float euclidian_6(feature_desc_6 x1, feature_desc_6 x2)
{

	float ab = abs(x1.s1 - x2.s1) / (x1.s1 < x2.s1 ? x1.s1 : x2.s1) +
		abs(x1.s2 - x2.s2) / (x1.s2 < x2.s2 ? x1.s2 : x2.s2) +
		abs(x1.s3 - x2.s3) / (x1.s3 < x2.s3 ? x1.s3 : x2.s3) +
		abs(x1.s4 - x2.s4) / (x1.s4 < x2.s4 ? x1.s4 : x2.s4) +
		abs(x1.s5 - x2.s5) / (x1.s5 < x2.s5 ? x1.s5 : x2.s5) +
		abs(x1.s6 - x2.s6) / (x1.s6 < x2.s6 ? x1.s6 : x2.s6) +
		abs(x1.s7 - x2.s7) / (x1.s7 < x2.s7 ? x1.s7 : x2.s7) +
		abs(x1.s8 - x2.s8) / (x1.s8 < x2.s8 ? x1.s8 : x2.s8) +
		abs(x1.s9 - x2.s9) / (x1.s9 < x2.s9 ? x1.s9 : x2.s9) +
		abs(x1.s10 - x2.s10) / (x1.s10 < x2.s10 ? x1.s10 : x2.s10)+
		abs(x1.s11 - x2.s11) / (x1.s11 < x2.s11 ? x1.s11 : x2.s11) +
		abs(x1.s12 - x2.s12) / (x1.s12 < x2.s12 ? x1.s12 : x2.s12) +
		abs(x1.s13 - x2.s13) / (x1.s13 < x2.s13 ? x1.s13 : x2.s13) +
		abs(x1.s14 - x2.s14) / (x1.s14 < x2.s14 ? x1.s14 : x2.s14) +
		abs(x1.s15 - x2.s15) / (x1.s15 < x2.s15 ? x1.s15 : x2.s15);

	return  ab;
}


int main(int argc, char** argv)
{
	ofstream fout("match_mod_tar.txt");
	pcl::console::TicToc time;
	time.tic();

	//加载模板点云
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mod(new pcl::PointCloud<pcl::PointXYZI>);
	if (pcl::io::loadPCDFile<pcl::PointXYZI>("mod.pcd", *cloud_mod) == -1) //* load the file
	{
		PCL_ERROR("Couldn't read file model file  \n");
		return (-1);
	}

	std::cout << "Loaded " << cloud_mod->points.size() << " data points  for model！" << std::endl;
	//加载模板点云

	//加载目标点云
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_tar(new pcl::PointCloud<pcl::PointXYZI>);
	if (pcl::io::loadPCDFile<pcl::PointXYZI>("tar.pcd", *cloud_tar) == -1) //* load the file
	{
		PCL_ERROR("Couldn't read file target file \n");
		return (-1);
	}

	std::cout << "Loaded " << cloud_tar->points.size() << " data points  for target ！" << std::endl;
	//加载目标点云

	fout << "load  mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
	//time.tic();

	//过滤点云
	time.tic();
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_filt_mod(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_filt_tar(new pcl::PointCloud<pcl::PointXYZI>);
	filter_intensity(cloud_mod, cloud_filt_mod, curv_thresh);
	filter_intensity(cloud_tar, cloud_filt_tar, curv_thresh);

	fout << "filter mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;

	pcl::io::savePCDFileASCII("cloud_filt_mod.pcd", *cloud_filt_mod);
	std::cout << "write point cloud into cloud_filt_mod!" << endl;
	pcl::io::savePCDFileASCII("cloud_filt_tar.pcd", *cloud_filt_tar);
	std::cout << "write point cloud into cloud_filt_tar!" << endl;

	time.tic();
	//过滤点云
	//生成特征点云团
	pcl::PointCloud<pcl::PointXYZI>::Ptr cluster_cloud_mod(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cluster_cloud_tar(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mean_points_mod(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mean_points_tar(new pcl::PointCloud<pcl::PointXYZI>);

	cluster_mean(cloud_filt_mod, &cluster_cloud_mod, &cloud_mean_points_mod);
	cluster_mean(cloud_filt_tar, &cluster_cloud_tar, &cloud_mean_points_tar);


	fout << "generate clusters and mean points of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;

	pcl::io::savePCDFileASCII("cluster_cloud_mod.pcd", *cluster_cloud_mod);
	std::cout << "write point cloud into cluster_cloud_mod!" << endl;
	pcl::io::savePCDFileASCII("cloud_mean_points_mod.pcd", *cloud_mean_points_mod);
	std::cout << "write point cloud into cloud_mean_points_mod!" << endl;
	pcl::io::savePCDFileASCII("cluster_cloud_tar.pcd", *cluster_cloud_tar);
	std::cout << "write point cloud into cluster_cloud_tar!" << endl;
	pcl::io::savePCDFileASCII("cloud_mean_points_tar.pcd", *cloud_mean_points_tar);
	std::cout << "write point cloud into cloud_mean_points_tar!" << endl;

	time.tic();
	//计算特征描述子

	std::vector<point_pair> match_mod_tar;
	std::vector<point_pair> match_tar_mod;


	
	if (desc_type==4) {
		std::vector<feature_desc_4> fea_vec_mod;
		std::vector<feature_desc_4> fea_vec_tar;

		feature_discriptor_4(fea_vec_mod, cloud_mean_points_mod);
		feature_discriptor_4(fea_vec_tar, cloud_mean_points_tar);

		for (std::vector<feature_desc_4>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;

		}

		for (std::vector<feature_desc_4>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " tar: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;

		}

		cout << "load files  1 " << endl;


		fout << "generate descriptors of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();

		//计算特征描述子

		//计算相似度

		// 计算 tar - mod最佳匹配 


		cout << "load files 232323 " << endl;
		for (std::vector<feature_desc_4>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			float min_rel = 50.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc_4>::const_iterator ir = fea_vec_mod.begin(); ir != fea_vec_mod.end(); ir++)
			{
				float rel_temp = euclidian_4(*it, *ir);

				if (rel_temp < 100)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << "  " << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_tar_mod.push_back(temp_pair);
		}
		cout << "load files 323232 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		//按误差由小到大 排序
		sort(match_tar_mod.begin(), match_tar_mod.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		cout << "load files 2 " << endl;
		for (std::vector<feature_desc_4>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			float min_rel = 100.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc_4>::const_iterator ir = fea_vec_tar.begin(); ir != fea_vec_tar.end(); ir++)
			{
				float rel_temp = euclidian_4(*it, *ir);

				if (rel_temp < 100)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_mod_tar.push_back(temp_pair);
		}
		cout << "load files 3 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}

		//按误差由小到大 排序
		sort(match_mod_tar.begin(), match_mod_tar.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}
		//todo

	std::vector<point_pair> mutual_match_mod_tar;
	for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); it++)
	{
	for (std::vector<point_pair>::const_iterator ir = match_mod_tar.begin(); ir != match_mod_tar.end(); ir++)
	{
	if (it->idx_2 == ir->idx_1&&it->idx_1 == ir->idx_2)
	{
	point_pair temp_pair;
	temp_pair.idx_1 = it->idx_2;
	temp_pair.idx_2 = it->idx_1;
	temp_pair.d = it->d;
	//cout << "bingo!  " << endl;
	fout << "bingo!  " << endl;
	//cout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;
	fout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;

	mutual_match_mod_tar.push_back(temp_pair);
	}
	}

	}
	sort(mutual_match_mod_tar.begin(), mutual_match_mod_tar.end(), LessSort);

	pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures(new pcl::PointCloud<pcl::PrincipalCurvatures>());

	//mod-tar ==tar-mod

	pcl::TransformationFromCorrespondences transformationFromCorr;

	for (std::vector<point_pair>::const_iterator it = mutual_match_mod_tar.begin(); it != mutual_match_mod_tar.begin() + reg_size; ++it)
	{
	Eigen::Vector3f from(cloud_mean_points_tar->points[it->idx_2].x,
	cloud_mean_points_tar->points[it->idx_2].y,
	cloud_mean_points_tar->points[it->idx_2].z);
	Eigen::Vector3f  to(cloud_mean_points_mod->points[it->idx_1].x,
	cloud_mean_points_mod->points[it->idx_1].y,
	cloud_mean_points_mod->points[it->idx_1].z);
	transformationFromCorr.add(from, to, 1.0);//all the same weight

	}


	fout << "generate match pairs of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
	time.tic();


	Eigen::Matrix4f transformationCorrespondence;
	transformationCorrespondence = transformationFromCorr.getTransformation().matrix();
	std::cout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;
	fout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;


	fout << "compute transform matrix of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;

	fout << "parameter ---- search_range: " << search_range << "MIU: " << MIU << "; SIGMA: " << SIGMA << "; min_cluster: " << min_cluster << "; curv_thresh: " << curv_thresh << "; reg_size: " << reg_size << endl;

	fout.close();

	}
	
	else if (desc_type == 3) {
		std::vector<feature_desc> fea_vec_mod;
		std::vector<feature_desc> fea_vec_tar;

		feature_discriptor(fea_vec_mod, cloud_mean_points_mod);
		feature_discriptor(fea_vec_tar, cloud_mean_points_tar);

		for (std::vector<feature_desc>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;

		}

		for (std::vector<feature_desc>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " tar: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;

		}

		cout << "load files  1 " << endl;


		fout << "generate descriptors of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();
		cout << "load files 232323 " << endl;
		for (std::vector<feature_desc>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			float min_rel = 50.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc>::const_iterator ir = fea_vec_mod.begin(); ir != fea_vec_mod.end(); ir++)
			{
				float rel_temp = euclidian(*it, *ir);

				if (rel_temp < 100)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << "  " << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_tar_mod.push_back(temp_pair);
		}
		cout << "load files 323232 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		//按误差由小到大 排序
		sort(match_tar_mod.begin(), match_tar_mod.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		cout << "load files 2 " << endl;
		for (std::vector<feature_desc>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			float min_rel = 100.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc>::const_iterator ir = fea_vec_tar.begin(); ir != fea_vec_tar.end(); ir++)
			{
				float rel_temp = euclidian(*it, *ir);

				if (rel_temp < 100)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_mod_tar.push_back(temp_pair);
		}
		cout << "load files 3 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}


		//按误差由小到大 排序
		sort(match_mod_tar.begin(), match_mod_tar.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}

		//fout.close();

		// mod-tar ==tar-mod


		std::vector<point_pair> mutual_match_mod_tar;
		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); it++)
		{
			for (std::vector<point_pair>::const_iterator ir = match_mod_tar.begin(); ir != match_mod_tar.end(); ir++)
			{
				if (it->idx_2 == ir->idx_1&&it->idx_1 == ir->idx_2)
				{
					point_pair temp_pair;
					temp_pair.idx_1 = it->idx_2;
					temp_pair.idx_2 = it->idx_1;
					temp_pair.d = it->d;
					//cout << "bingo!  " << endl;
					fout << "bingo!  " << endl;
					//cout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;
					fout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;

					mutual_match_mod_tar.push_back(temp_pair);
				}
			}

		}
		sort(mutual_match_mod_tar.begin(), mutual_match_mod_tar.end(), LessSort);

		pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures(new pcl::PointCloud<pcl::PrincipalCurvatures>());

		pcl::TransformationFromCorrespondences transformationFromCorr;

		for (std::vector<point_pair>::const_iterator it = mutual_match_mod_tar.begin(); it != mutual_match_mod_tar.begin() + reg_size; ++it)
		{
			Eigen::Vector3f from(cloud_mean_points_tar->points[it->idx_2].x,
				cloud_mean_points_tar->points[it->idx_2].y,
				cloud_mean_points_tar->points[it->idx_2].z);
			Eigen::Vector3f  to(cloud_mean_points_mod->points[it->idx_1].x,
				cloud_mean_points_mod->points[it->idx_1].y,
				cloud_mean_points_mod->points[it->idx_1].z);
			transformationFromCorr.add(from, to, 1.0);//all the same weight 

		}

		fout << "generate match pairs of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();


		Eigen::Matrix4f transformationCorrespondence;
		transformationCorrespondence = transformationFromCorr.getTransformation().matrix();
		std::cout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;
		fout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;


		fout << "compute transform matrix of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;

		fout << "parameter ---- search_range: " << search_range << "MIU: " << MIU << "; SIGMA: " << SIGMA << "; min_cluster: " << min_cluster << "; curv_thresh: " << curv_thresh << "; reg_size: " << reg_size << endl;

		fout.close();
	}
	
	else if (desc_type == 5) {
		std::vector<feature_desc_5> fea_vec_mod;
		std::vector<feature_desc_5> fea_vec_tar;

		feature_discriptor_5(fea_vec_mod, cloud_mean_points_mod);
		feature_discriptor_5(fea_vec_tar, cloud_mean_points_tar);

		for (std::vector<feature_desc_5>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;

		}

		for (std::vector<feature_desc_5>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " tar: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;

		}

		cout << "load files  1 " << endl;

		fout << "generate descriptors of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();

		cout << "load files 232323 " << endl;
		for (std::vector<feature_desc_5>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			float min_rel = 50.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc_5>::const_iterator ir = fea_vec_mod.begin(); ir != fea_vec_mod.end(); ir++)
			{
				float rel_temp = euclidian_5(*it, *ir);

				if (rel_temp < 10)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << "  " << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_tar_mod.push_back(temp_pair);
		}
		cout << "load files 323232 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		//按误差由小到大 排序
		sort(match_tar_mod.begin(), match_tar_mod.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		cout << "load files 2 " << endl;
		for (std::vector<feature_desc_5>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			float min_rel = 100.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc_5>::const_iterator ir = fea_vec_tar.begin(); ir != fea_vec_tar.end(); ir++)
			{
				float rel_temp = euclidian_5(*it, *ir);

				if (rel_temp < 10)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_mod_tar.push_back(temp_pair);
		}
		cout << "load files 3 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}


		//按误差由小到大 排序
		sort(match_mod_tar.begin(), match_mod_tar.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}
		//todo

		std::vector<point_pair> mutual_match_mod_tar;
		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); it++)
		{
			for (std::vector<point_pair>::const_iterator ir = match_mod_tar.begin(); ir != match_mod_tar.end(); ir++)
			{
				if (it->idx_2 == ir->idx_1&&it->idx_1 == ir->idx_2)
				{
					point_pair temp_pair;
					temp_pair.idx_1 = it->idx_2;
					temp_pair.idx_2 = it->idx_1;
					temp_pair.d = it->d;
					//cout << "bingo!  " << endl;
					fout << "bingo!  " << endl;
					//cout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;
					fout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;

					mutual_match_mod_tar.push_back(temp_pair);
				}
			}

		}
		sort(mutual_match_mod_tar.begin(), mutual_match_mod_tar.end(), LessSort);

		pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures(new pcl::PointCloud<pcl::PrincipalCurvatures>());

		//mod-tar ==tar-mod

		pcl::TransformationFromCorrespondences transformationFromCorr;

		for (std::vector<point_pair>::const_iterator it = mutual_match_mod_tar.begin(); it != mutual_match_mod_tar.begin() + reg_size; ++it)
		{
			Eigen::Vector3f from(cloud_mean_points_tar->points[it->idx_2].x,
				cloud_mean_points_tar->points[it->idx_2].y,
				cloud_mean_points_tar->points[it->idx_2].z);
			Eigen::Vector3f  to(cloud_mean_points_mod->points[it->idx_1].x,
				cloud_mean_points_mod->points[it->idx_1].y,
				cloud_mean_points_mod->points[it->idx_1].z);
			transformationFromCorr.add(from, to, 1.0);//all the same weight

		}
		fout << "generate match pairs of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();

		Eigen::Matrix4f transformationCorrespondence;
		transformationCorrespondence = transformationFromCorr.getTransformation().matrix();
		std::cout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;
		fout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;

		fout << "compute transform matrix of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		fout << "parameter ---- search_range: " << search_range << "MIU: " << MIU << "; SIGMA: " << SIGMA << "; min_cluster: " << min_cluster << "; curv_thresh: " << curv_thresh << "; reg_size: " << reg_size << " ; desc_type: "<< desc_type << endl;

		fout.close();

	}

	else if (desc_type == 6) {
		std::vector<feature_desc_6> fea_vec_mod;
		std::vector<feature_desc_6> fea_vec_tar;

		feature_discriptor_6(fea_vec_mod, cloud_mean_points_mod);
		feature_discriptor_6(fea_vec_tar, cloud_mean_points_tar);

		for (std::vector<feature_desc_6>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
		}

		for (std::vector<feature_desc_6>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			//cout << " mod: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
			fout << " tar: " << it->idx << "  " << it->s1 << "  " << it->s2 << "  " << it->s3 << "  " << it->i1 << "  " << it->i2 << "  " << it->i3 << endl;
		}

		cout << "load files  1 " << endl;


		fout << "generate descriptors of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();

		cout << "load files 232323 " << endl;
		for (std::vector<feature_desc_6>::const_iterator it = fea_vec_tar.begin(); it != fea_vec_tar.end(); it++)
		{
			float min_rel = 50.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc_6>::const_iterator ir = fea_vec_mod.begin(); ir != fea_vec_mod.end(); ir++)
			{
				float rel_temp = euclidian_6(*it, *ir);

				if (rel_temp < 10)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " tar: " << it->idx << " " << cloud_mean_points_tar->points[it->idx].x << " " << cloud_mean_points_tar->points[it->idx].y << " " << cloud_mean_points_tar->points[it->idx].z << "  " << cloud_mean_points_tar->points[it->idx].intensity << " " << " " << " mod: " << ir->idx << " " << cloud_mean_points_mod->points[ir->idx].x << " " << cloud_mean_points_mod->points[ir->idx].y << " " << cloud_mean_points_mod->points[ir->idx].z << " " << cloud_mean_points_mod->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_tar_mod.push_back(temp_pair);
		}
		cout << "load files 323232 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		//按误差由小到大 排序
		sort(match_tar_mod.begin(), match_tar_mod.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); ++it)
		{
			//	cout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
			fout << " tar: " << it->idx_1 << " mod: " << it->idx_2 << " res: " << it->d << endl;
		}

		cout << "load files 2 " << endl;
		for (std::vector<feature_desc_6>::const_iterator it = fea_vec_mod.begin(); it != fea_vec_mod.end(); it++)
		{
			float min_rel = 100.0;//最小误差
			int  arch = -1;//最小误差匹配位置

			for (std::vector<feature_desc_6>::const_iterator ir = fea_vec_tar.begin(); ir != fea_vec_tar.end(); ir++)
			{
				float rel_temp = euclidian_6(*it, *ir);

				if (rel_temp < 10)
				{
					if (min_rel > rel_temp)
					{
						min_rel = rel_temp;
						arch = ir->idx;
					}
					//	cout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;
					fout << " mod: " << it->idx << " " << cloud_mean_points_mod->points[it->idx].x << " " << cloud_mean_points_mod->points[it->idx].y << " " << cloud_mean_points_mod->points[it->idx].z << cloud_mean_points_mod->points[it->idx].intensity << " " << " " << " tar: " << ir->idx << " " << cloud_mean_points_tar->points[ir->idx].x << " " << cloud_mean_points_tar->points[ir->idx].y << " " << cloud_mean_points_tar->points[ir->idx].z << " " << cloud_mean_points_tar->points[ir->idx].intensity << " " << " res: " << rel_temp << endl;

					fout << " mod: " << it->idx << " tar: " << ir->idx << " res: " << rel_temp << endl;

				}
			}
			point_pair temp_pair;
			temp_pair.idx_1 = it->idx;
			temp_pair.idx_2 = arch;
			temp_pair.d = min_rel;

			match_mod_tar.push_back(temp_pair);
		}
		cout << "load files 3 " << endl;
		//	ofstream fout("match_mod_tar.txt");

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}

		//按误差由小到大 排序
		sort(match_mod_tar.begin(), match_mod_tar.end(), LessSort);//

		for (std::vector<point_pair>::const_iterator it = match_mod_tar.begin(); it != match_mod_tar.end(); ++it)
		{
			//	cout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
			fout << " mod: " << it->idx_1 << " tar: " << it->idx_2 << " res: " << it->d << endl;
		}
		//todo

		std::vector<point_pair> mutual_match_mod_tar;
		for (std::vector<point_pair>::const_iterator it = match_tar_mod.begin(); it != match_tar_mod.end(); it++)
		{
			for (std::vector<point_pair>::const_iterator ir = match_mod_tar.begin(); ir != match_mod_tar.end(); ir++)
			{
				if (it->idx_2 == ir->idx_1&&it->idx_1 == ir->idx_2)
				{
					point_pair temp_pair;
					temp_pair.idx_1 = it->idx_2;
					temp_pair.idx_2 = it->idx_1;
					temp_pair.d = it->d;
					//cout << "bingo!  " << endl;
					fout << "bingo!  " << endl;
					//cout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;
					fout << " mod: " << it->idx_2 << " tar: " << it->idx_1 << " res: " << it->d << endl;

					mutual_match_mod_tar.push_back(temp_pair);
				}
			}

		}
		sort(mutual_match_mod_tar.begin(), mutual_match_mod_tar.end(), LessSort);

		pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures(new pcl::PointCloud<pcl::PrincipalCurvatures>());

		//mod-tar ==tar-mod

		pcl::TransformationFromCorrespondences transformationFromCorr;

		for (std::vector<point_pair>::const_iterator it = mutual_match_mod_tar.begin(); it != mutual_match_mod_tar.begin() + reg_size; ++it)
		{
			Eigen::Vector3f from(cloud_mean_points_tar->points[it->idx_2].x,
				cloud_mean_points_tar->points[it->idx_2].y,
				cloud_mean_points_tar->points[it->idx_2].z);
			Eigen::Vector3f  to(cloud_mean_points_mod->points[it->idx_1].x,
				cloud_mean_points_mod->points[it->idx_1].y,
				cloud_mean_points_mod->points[it->idx_1].z);
			transformationFromCorr.add(from, to, 1.0);//all the same weight
		}

		fout << "generate match pairs of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		time.tic();

		Eigen::Matrix4f transformationCorrespondence;
		transformationCorrespondence = transformationFromCorr.getTransformation().matrix();
		std::cout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;
		fout << "\ntransformation from corresponding points is \n" << transformationCorrespondence << std::endl;
		fout << "compute transform matrix of mod and tar files consume time: " << time.toc() / 1000 << "s" << endl;
		fout << "parameter ---- search_range: " << search_range << "MIU: " << MIU << "; SIGMA: " << SIGMA << "; min_cluster: " << min_cluster << "; curv_thresh: " << curv_thresh << "; reg_size: " << reg_size << " ; desc_type: " << desc_type << endl;
		fout.close();

	}
	
}