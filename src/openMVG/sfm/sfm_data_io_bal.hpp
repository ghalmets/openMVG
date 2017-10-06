
// Copyright (c) 2015 Pierre MOULON.

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

void RotationMatrixToQuaternion(
		const Eigen::Matrix3d& R,
		double* quaternion) {
	const double trace = R(0, 0) + R(1, 1) + R(2, 2);
	if (trace >= 0.0) {
		double t = sqrt(trace + double(1.0));
		quaternion[0] = double(0.5) * t;
		t = double(0.5) / t;
		quaternion[1] = (R(2, 1) - R(1, 2)) * t;
		quaternion[2] = (R(0, 2) - R(2, 0)) * t;
		quaternion[3] = (R(1, 0) - R(0, 1)) * t;
	} else {
		int i = 0;
		if (R(1, 1) > R(0, 0)) {
			i = 1;
		}

		if (R(2, 2) > R(i, i)) {
			i = 2;
		}

		const int j = (i + 1) % 3;
		const int k = (j + 1) % 3;
		double t = sqrt(R(i, i) - R(j, j) - R(k, k) + double(1.0));
		quaternion[i + 1] = double(0.5) * t;
		t = double(0.5) / t;
		quaternion[0] = (R(k, j) - R(j, k)) * t;
		quaternion[j + 1] = (R(j, i) + R(i, j)) * t;
		quaternion[k + 1] = (R(k, i) + R(i, k)) * t;
	}
}

inline void QuaternionToAngleAxis(const double* quaternion, double* angle_axis) {
	const double& q1 = quaternion[1];
	const double& q2 = quaternion[2];
	const double& q3 = quaternion[3];
	const double sin_squared_theta = q1 * q1 + q2 * q2 + q3 * q3;

	// For quaternions representing non-zero rotation, the conversion
	// is numerically stable.
	if (sin_squared_theta > double(0.0)) {
		const double sin_theta = sqrt(sin_squared_theta);
		const double& cos_theta = quaternion[0];

		// If cos_theta is negative, theta is greater than pi/2, which
		// means that angle for the angle_axis vector which is 2 * theta
		// would be greater than pi.
		//
		// While this will result in the correct rotation, it does not
		// result in a normalized angle-axis vector.
		//
		// In that case we observe that 2 * theta ~ 2 * theta - 2 * pi,
		// which is equivalent saying
		//
		//   theta - pi = atan(sin(theta - pi), cos(theta - pi))
		//              = atan(-sin(theta), -cos(theta))
		//
		const double two_theta =
				double(2.0) * ((cos_theta < 0.0)
							   ? atan2(-sin_theta, -cos_theta)
							   : atan2(sin_theta, cos_theta));
		const double k = two_theta / sin_theta;
		angle_axis[0] = q1 * k;
		angle_axis[1] = q2 * k;
		angle_axis[2] = q3 * k;
	} else {
		// For zero rotation, sqrt() will produce NaN in the derivative since
		// the argument is zero.  By approximating with a Taylor series,
		// and truncating at one term, the value and first derivatives will be
		// computed correctly when Jets are used.
		const double k(2.0);
		angle_axis[0] = q1 * k;
		angle_axis[1] = q2 * k;
		angle_axis[2] = q3 * k;
	}
}

void RotationMatrixToAngleAxis(
		Eigen::Matrix3d& R,
		double* angle_axis) {
	double quaternion[4];
	RotationMatrixToQuaternion(R, quaternion);
	QuaternionToAngleAxis(quaternion, angle_axis);
	return;
}


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

		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];
		time (&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer,80,"%Y%m%d%__H%M%S",timeinfo);
		std::string time(buffer);

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

			write	<< 0 << "\n" ; //f
			write	<<	0 << "\n" ; //k1
			write	<<	0 << "\n" ; //k2
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
			//OpenMV provides us with C (Camera position from Global COS, , but we need t, so that uv=K*(R|t)X;
						camera_poses[i][3]=-t[0];
						camera_poses[i][4]=-t[1];
						camera_poses[i][5]=-t[2];

//			camera_poses[i][3]=rot_inv[i][3];
//			camera_poses[i][4]=rot_inv[i][7];
//			camera_poses[i][5]=rot_inv[i][11];

			//			Eigen::Matrix3d rot2rodrigues;
			//			rot2rodrigues<<	rot_inv(0,0),rot_inv(0,1),rot_inv(0,2),
			//					rot_inv(1,0),rot_inv(1,1),rot_inv(1,2),
			//					rot_inv(2,0),rot_inv(2,1),rot_inv(2,2);

			Eigen::Matrix3d rot2rodrigues;
			rot2rodrigues<<camera_poses_[i][0],		camera_poses_[i][1],	camera_poses_[i][2],
					camera_poses_[i][4],	camera_poses_[i][5],	camera_poses_[i][6],
					camera_poses_[i][8],	camera_poses_[i][9],	camera_poses_[i][10];

//			rot2rodrigues<<	rot(0,0),rot(0,1),rot(0,2),
//					rot(1,0),rot(1,1),rot(1,2),
//					rot(2,0),rot(2,1),rot(2,2);

			//			Eigen::EigenSolver<Eigen::Matrix3d > solver(rot2rodrigues,true);
			//			Eigen::Vector3d angle_axis;

			//			double angle_axis[3];

			//			RotationMatrixToAngleAxis(rot2rodrigues,angle_axis);

			//			camera_poses[i][0]=angle_axis[0];
			//			camera_poses[i][1]=angle_axis[1];
			//			camera_poses[i][2]=angle_axis[2];

			Eigen::AngleAxis<double >rodrig;
			rodrig.fromRotationMatrix(rot2rodrigues);

			//Convert 4param Rodrig to 3param  Rodrigues

			camera_poses[i][0]=rodrig.axis()[0]*rodrig.angle();
			camera_poses[i][1]=rodrig.axis()[1]*rodrig.angle();
			camera_poses[i][2]=rodrig.axis()[2]*rodrig.angle();





			//			Eigen::AngleAxis<double >rodrig;
			//			rodrig.fromRotationMatrix(rot2rodrigues);

			//			//Convert Quaternion to 3param  Rodrigues

			//			camera_poses[i*6+0]=rodrig.axis()[0]*rodrig.angle();
			//			camera_poses[i*6+1]=rodrig.axis()[1]*rodrig.angle();
			//			camera_poses[i*6+2]=rodrig.axis()[2]*rodrig.angle();


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

/// Save SfM_Data in an ASCII BAF (Bundle Adjustment File).
// --Header
// #Intrinsics
// #Poses
// #Landmarks
// --Data
// Intrinsic parameters [foc ppx ppy, ...]
// Poses [angle axis, camera center]
// Landmarks [X Y Z #observations id_intrinsic id_pose x y ...]
//--
//- Export also a _imgList.txt file with View filename and id_intrinsic & id_pose.
// filename id_intrinsic id_pose
// The ids allow to establish a link between 3D point observations & the corresponding views
//--
// Export missing poses as Identity pose to keep tracking of the original id_pose indexes
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

	//double* camera_poses=new double[num_views*12];
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

			if(view_num >=num_views ){
				std::cout<<"VIEW_KEY EXCEEDS NUM_VIEWS!!!"<<std::endl;

			}

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

				//Row Major
				//				camera_poses[view_num][0]=rotation[0];
				//				camera_poses[view_num][1]=rotation[1];
				//				camera_poses[view_num][2]=rotation[2];

				//				camera_poses[view_num][4]=rotation[3];
				//				camera_poses[view_num][5]=rotation[4];
				//				camera_poses[view_num][6]=rotation[5];

				//				camera_poses[view_num][8]=rotation[6];
				//				camera_poses[view_num][9]=rotation[7];
				//				camera_poses[view_num][10]=rotation[8];

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

		std::cout<<"Exporting Landmarks"<<std::endl;
		const Landmarks & landmarks = sfm_data.GetLandmarks();

		//unsigned int length=3*num_landmarks;
		//double* keypoints=new double[length];
		std::vector <std::vector <double> >keypoints(num_landmarks, std::vector<double>(3));


		std::cout<<"keypoints generated on heap"<<std::endl;
		//Export Features
		unsigned int i=0;
		for (const auto & iterLandmarks : landmarks )
		{
			// Export visibility information
			// X Y Z #observations id_cam id_pose x y ...
			const double * X = iterLandmarks.second.X.data();
			unsigned int feature_key=iterLandmarks.first;

			if(feature_key>=num_landmarks){
				std::cout<<"FEATURE KEY EXCEEDS ARRAY SIZE!!! "<<feature_key<<i<<std::endl;
			}
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

				std::cout<<"[feature_id,view_id, [u,v], pose_id, [pose], image_path]: "<<i<<" , "<<view_key<<" , "<<"["<<uv[0]<<" ,  "<<uv[1]<<"]"<<", "<< v->id_pose;

				const double * rotation = poses.at(v->id_pose).rotation().data();
				const double * center = poses.at(v->id_pose).center().data();

				std::cout<<", [ "<<center[0]<<" , "
						<<center[1]<<" , "
					   <<center[2]<<" , "
					  <<rotation[0]<<" , "
					 <<rotation[1]<<" , "
					<<rotation[2]<<" , "
				   <<rotation[3]<<" , "
				  <<rotation[4]<<" , "
				 <<rotation[5]<< " ,"
				<<rotation[6]<<" , "
			   <<rotation[7]<<" , "
			  <<rotation[8]<<" ]"
				<< v->s_Img_path<<std::endl;



				//				const View * view = iterOb.first;
				//				int view_num = view->id_view;

				//				double uv[2];
				//				uv[0]= iterOb.second.x(0);
				//				uv[1]= iterOb.second.x(1);
				//				observations.addObservation(view_num,i, uv);
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

	// Export depth data filenames & ids as an imgList.txt file
	{

		const char* path = filename.c_str();
		boost::filesystem::path file_path(path);
		boost::filesystem::path dir_d(file_path.parent_path().string()+"/BAL/depth/");
		if(boost::filesystem::create_directory(dir_d))
		{
			std::cout<< "Depth Directory Created: "<<dir_d.string()<<std::endl;
		}

		boost::filesystem::path file_path_rgb(path);
		boost::filesystem::path dir_rgb(file_path_rgb.parent_path().string()+"/BAL/rgb/");
		if(boost::filesystem::create_directory(dir_rgb))
		{
			std::cout<< "RGB Directory Created: "<<dir_rgb.string()<<std::endl;
		}


		boost::filesystem::path dir_img_ori=file_path.parent_path();

		std::cout<<"Reading depth/rgb data from"<<dir_img_ori.string()<<"/../../[depth|rgb]/"<<std::endl;


		for(size_t i=0;i<camera_path.size();i++)	{


			std::stringstream ss_depth_in;
			ss_depth_in<<dir_img_ori.string()<<"/../../depth/"<<camera_path[i];
			std::ifstream  src_d(ss_depth_in.str(), std::ios::binary);

			std::stringstream ss_rgb_in;
			ss_rgb_in<<dir_img_ori.string()<<"/../../rgb/"<<camera_path[i];
			std::ifstream  src_rgb(ss_rgb_in.str(), std::ios::binary);

			//			std::stringstream ss_d_out;
			//			ss_d_out<<dir_d.string()<<iterV.second->id_view<<".png";
			//			std::ofstream  dst_d(ss_d_out.str(),   std::ios::binary);

			//			std::stringstream ss_rgb_out;
			//			ss_rgb_out<<dir_rgb.string()<<iterV.second->id_view<<".png";
			//			std::ofstream  dst_rgb(ss_rgb_out.str(),   std::ios::binary);

			std::stringstream ss_d_out;
			ss_d_out<<dir_d.string()<<i<<".png";
			std::ofstream  dst_d(ss_d_out.str(),   std::ios::binary);

			std::stringstream ss_rgb_out;
			ss_rgb_out<<dir_rgb.string()<<i<<".png";
			std::ofstream  dst_rgb(ss_rgb_out.str(),   std::ios::binary);

			dst_d << src_d.rdbuf();
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



// This algorithm comes from "Quaternion Calculus and Fast Animation",
// Ken Shoemake, 1987 SIGGRAPH course notes

#endif // OPENMVG_SFM_SFM_DATA_IO_BAL_HPP
