
// Copyright (c) 2015 Pierre MOULON.
// Copyright (c) 2017 Georg Halmetschlager-Funek, gh@acin.tuwien.ac.at
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_IO_BAL_HPP
#define OPENMVG_SFM_SFM_DATA_IO_BAL_HPP

#include "openMVG/sfm/sfm_data_io.hpp"

#include <fstream>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>

#include <boost/filesystem.hpp>

#define EXPORT_DEPTH

namespace openMVG {
namespace sfm {


class generateBALObservations{

public:

	~generateBALObservations(){

	}

	void addObservation(int camera_index, int key_index, double* uv){
		BALObservation add;
		add.camera_index=camera_index;
		add.key_index=key_index;
		add.uv[0]=uv[0];
		add.uv[1]=uv[1];
		observations.push_back(add);

	}
	void writeToFile(std::string path){
		std::fstream write;
		double num_observations_;

		write.open(path+"/bal.data", std::fstream::out);

		num_observations_= double(observations.size()) ;

		write << num_cameras<<" "<< num_keypoints <<" "<< num_observations_<<"\n";

		for(size_t i=0; i<observations.size(); i++){ //Write observations

			write	<<	observations[i].camera_index <<" "<< observations[i].key_index <<" "<< observations[i].uv[0] <<" "<< observations[i].uv[1] <<"\n";
		}

		for(size_t i=0; i<num_cameras; i++){ //Write Camera Poses and

			write	<<	camera_poses[i][0] << "\n" ;
			write	<<	camera_poses[i][1] << "\n" ;
			write	<<	camera_poses[i][2] << "\n" ;
			write	<<	camera_poses[i][3] << "\n" ;
			write	<<	camera_poses[i][4] << "\n" ;
			write	<<	camera_poses[i][5] << "\n" ;

			write	<< 0.0 << "\n" ; //f, not exported
			write	<<	0.0 << "\n" ; //k1, not exported
			write	<<	0.0 << "\n" ; //k2, not exported
		}
		for(size_t i=0; i<num_keypoints; i++){
			for(size_t j=0; j<3; j++){
				write	<<	key_points[i][j] << "\n" ;
			}
		}
		write.close();

	}

	void setNumCameras(unsigned int numcameras_){num_cameras=numcameras_;}
	void setNumKeypoints(unsigned int num_keypoints_){num_keypoints=num_keypoints_;}
	void setCameraPoses(std::vector< std::vector<double> > camera_poses_)
	{
		camera_poses.resize(num_cameras);
		for(size_t i=0;i<camera_poses.size();i++){
			camera_poses[i].resize(6);
		}

		for(size_t i=0; i<num_cameras;i++){


			Eigen::Matrix3d rot;
			rot<<camera_poses_[i][0],		camera_poses_[i][1],	camera_poses_[i][2],
					camera_poses_[i][4],	camera_poses_[i][5],	camera_poses_[i][6],
					camera_poses_[i][8],	camera_poses_[i][9],	camera_poses_[i][10];


			Eigen::Vector3d c;
			Eigen::Vector3d t;

			c[0]=	camera_poses_[i][3];
			c[1]=	camera_poses_[i][7];
			c[2]=	camera_poses_[i][11];

			t=rot*c;
			//OpenMVG provides C (Camera position from Global COS, but we need t, so that uv=K*(R|t)X;
			camera_poses[i][3]=-t[0];
			camera_poses[i][4]=-t[1];
			camera_poses[i][5]=-t[2];


			Eigen::Matrix3d rot2rodrigues;
			rot2rodrigues<<camera_poses_[i][0],		camera_poses_[i][1],	camera_poses_[i][2],
					camera_poses_[i][4],	camera_poses_[i][5],	camera_poses_[i][6],
					camera_poses_[i][8],	camera_poses_[i][9],	camera_poses_[i][10];

			Eigen::AngleAxis<double >rodrig;
			rodrig.fromRotationMatrix(rot2rodrigues);

			//Convert 4param Rodrig to 3param  Rodrigues
			camera_poses[i][0]=rodrig.axis()[0]*rodrig.angle();
			camera_poses[i][1]=rodrig.axis()[1]*rodrig.angle();
			camera_poses[i][2]=rodrig.axis()[2]*rodrig.angle();

		}

	}
	void setKeypointPosition(std::vector< std::vector< double> > keypts){

		key_points=keypts;
	}

private:
	struct BALObservation{
		unsigned int camera_index;
		unsigned int key_index;
		double uv[2];
	};
	std::vector< std::vector< double> > camera_poses;
	std::vector< std::vector< double> > key_points;
	unsigned int num_cameras;
	unsigned int num_keypoints;
	std::vector<BALObservation> observations;

};

/// Save SfM_Data in an ASCII BAL (Bundle Adjustment in the Large). And exports RGB (& depth files) to camera_index.png
//format: http://grail.cs.washington.edu/projects/bal/
//exports in addition all rgb and depth data ordered from 0...1


inline bool Save_BAL(
		const SfM_Data & sfm_data,
		const std::string & filename,
		ESfM_Data flags_part)
{
	std::ofstream stream(filename.c_str());
	if (!stream.is_open())
		return false;

	unsigned int num_landmarks=sfm_data.GetLandmarks().size();
	unsigned int num_views=sfm_data.GetViews().size();

	std::cout<<"Generating BAL Observations..."<<std::endl;
	generateBALObservations observations;

	observations.setNumCameras(num_views);
	observations.setNumKeypoints(num_landmarks );

	std::vector <std::vector<double> > camera_poses(num_views, std::vector<double>(12));
	std::vector < std::string > camera_path(num_views);

	bool bOk = false;
	{

		std::cout<<"Exporting Camera Poses"<<std::endl;


		//Export Camera Poses

		const Poses & poses = sfm_data.GetPoses();
		unsigned int view_num=0;
		for ( const auto & iterV : sfm_data.GetViews() )
		{
			const View * view = iterV.second.get();

			view_num=iterV.first;

			camera_path[view_num]=view->s_Img_path;

			if (!sfm_data.IsPoseAndIntrinsicDefined(view))
			{
				const Mat3 R = Mat3::Identity();
				const double * rotation = R.data();

				camera_poses[view_num][0]=rotation[0];
				camera_poses[view_num][1]=rotation[1];
				camera_poses[view_num][2]=rotation[2];

				camera_poses[view_num][4]=rotation[3];
				camera_poses[view_num][5]=rotation[4];
				camera_poses[view_num][6]=rotation[5];

				camera_poses[view_num][8]=rotation[6];
				camera_poses[view_num][9]=rotation[7];
				camera_poses[view_num][10]=rotation[8];


				const Vec3 C = Vec3::Zero();
				const double * center = C.data();
				camera_poses[view_num][3]=center[0];
				camera_poses[view_num][7]=center[1];
				camera_poses[view_num][11]=center[2];

			}
			else
			{
				// [Rotation col major 3x3; camera center 3x1]
				const double * rotation = poses.at(view->id_pose).rotation().data();

				//COL MAJOR
				camera_poses[view_num][0]=rotation[0];
				camera_poses[view_num][1]=rotation[3];
				camera_poses[view_num][2]=rotation[6];

				camera_poses[view_num][4]=rotation[1];
				camera_poses[view_num][5]=rotation[4];
				camera_poses[view_num][6]=rotation[7];

				camera_poses[view_num][8]=rotation[2];
				camera_poses[view_num][9]=rotation[5];
				camera_poses[view_num][10]=rotation[8];


				const double * center = poses.at(view->id_pose).center().data();
				camera_poses[view_num][3]=center[0];
				camera_poses[view_num][7]=center[1];
				camera_poses[view_num][11]=center[2];

			}
		}

		std::cout<<"Exporting Landmarks.."<<std::endl;
		const Landmarks & landmarks = sfm_data.GetLandmarks();

		std::vector <std::vector <double> >keypoints(num_landmarks, std::vector<double>(3));

		//Export Features
		unsigned int i=0;
		for (const auto & iterLandmarks : landmarks )
		{
			// Export visibility information
			// X Y Z #observations id_cam id_pose x y ...
			const double * X = iterLandmarks.second.X.data();
			unsigned int feature_key=iterLandmarks.first;

			keypoints[i][0]=X[0];
			keypoints[i][1]=X[1];
			keypoints[i][2]=X[2];

			//Export Observation

			const Observations & obs = iterLandmarks.second.obs;

			for ( const auto & iterOb : obs )
			{

				const IndexT view_key = iterOb.first;


				const View * v = sfm_data.GetViews().at(view_key).get();

				double uv[2];
				uv[0]= iterOb.second.x(0);
				uv[1]= iterOb.second.x(1);

				observations.addObservation(view_key, i, uv);
			}
			i++;

		}
		std::cout<<"Setting up observation"<<std::endl;

		std::cout<<"Camera Pose"<<std::endl;
		observations.setCameraPoses(camera_poses);

		std::cout<<"Keypoint"<<std::endl;
		observations.setKeypointPosition(keypoints);

		std::cout<<"Write Data.."<<std::endl;
		std::cout<<"Creating Directory './BAL/'"<<std::endl;

		const char* path = filename.c_str();
		boost::filesystem::path file_path(path);
		boost::filesystem::path dir(file_path.parent_path().string()+"/BAL/");
		if(boost::filesystem::create_directory(dir))
		{
			std::cout<< "Directory Created: "<<dir.string()<<std::endl;
		}

		observations.writeToFile(dir.string());

		stream.flush();
		bOk = stream.good();
		stream.close();
	}

//Export Images
	{

		const char* path = filename.c_str();
		boost::filesystem::path file_path(path);
#ifdef EXPORT_DEPTH
		boost::filesystem::path dir_d(file_path.parent_path().string()+"/BAL/depth/");
		if(boost::filesystem::create_directory(dir_d))
		{
			std::cout<< "Depth Directory Created: "<<dir_d.string()<<std::endl;
		}
#endif

		boost::filesystem::path file_path_rgb(path);
		boost::filesystem::path dir_rgb(file_path_rgb.parent_path().string()+"/BAL/rgb/");
		if(boost::filesystem::create_directory(dir_rgb))
		{
			std::cout<< "RGB Directory Created: "<<dir_rgb.string()<<std::endl;
		}


		boost::filesystem::path dir_img_ori=file_path.parent_path();

		std::cout<<"Reading depth/rgb data from"<<dir_img_ori.string()<<"/../../[depth|rgb]/"<<std::endl;


		for(size_t i=0;i<camera_path.size();i++)	{

#ifdef EXPORT_DEPTH

			std::stringstream ss_depth_in;
			ss_depth_in<<dir_img_ori.string()<<"/../../depth/"<<camera_path[i];
			std::ifstream  src_d(ss_depth_in.str(), std::ios::binary);

			std::stringstream ss_d_out;
			ss_d_out<<dir_d.string()<<i<<".png";
			std::ofstream  dst_d(ss_d_out.str(),   std::ios::binary);

			dst_d << src_d.rdbuf();
#endif

			std::stringstream ss_rgb_in;
			ss_rgb_in<<dir_img_ori.string()<<"/../../rgb/"<<camera_path[i];
			std::ifstream  src_rgb(ss_rgb_in.str(), std::ios::binary);

			std::stringstream ss_rgb_out;
			ss_rgb_out<<dir_rgb.string()<<i<<".png";
			std::ofstream  dst_rgb(ss_rgb_out.str(),   std::ios::binary);

			dst_rgb << src_rgb.rdbuf();

		}
		stream.flush();
		bOk = stream.good();
		stream.close();
	}
	return bOk;
}

} // namespace sfm
} // namespace openMVG


#endif // OPENMVG_SFM_SFM_DATA_IO_BAL_HPP
