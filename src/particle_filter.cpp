/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "helper_functions.h"

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  num_particles = 1000;
  yaw_rate_limit = 1e-4;

  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta,std[2]);

  for (int i = 0; i < num_particles; i++){
    Particle new_particle;
    new_particle.id     = i;
    new_particle.x      = dist_x(gen);
    new_particle.y      = dist_y(gen);
    new_particle.theta  = dist_theta(gen);
    new_particle.weight = 1;
    particles.push_back(new_particle);
    weights.push_back(1.0);
  }
  
  is_initialized = true;
  return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0,std_pos[2]);

  double yawdt = yaw_rate * delta_t;
  double vdtheta;
  if (fabs(yaw_rate) > yaw_rate_limit){
    vdtheta = velocity / yaw_rate;
  }else{
    vdtheta = velocity * delta_t; //not really vdtheta, but micro-optimizing to not create another variable
  }
  for (int i = 0; i < num_particles; i++){
    if(fabs(yaw_rate) > yaw_rate_limit){
      particles[i].x += vdtheta * (sin(particles[i].theta + yawdt) - sin(particles[i].theta)) + dist_x(gen);
      particles[i].y += vdtheta * (cos(particles[i].theta) - cos(particles[i].theta + yawdt)) + dist_y(gen);
    }else{
      particles[i].x += vdtheta * cos(particles[i].theta) + dist_x(gen);
      particles[i].y += vdtheta * sin(particles[i].theta) + dist_y(gen);
    }
    particles[i].theta += yawdt + dist_theta(gen);
  }//for i over num_particles
  return;
}

vector<LandmarkObs> ParticleFilter::find_close_landmarks(double sensor_range, Particle &pcle, const Map &map_landmarks){
  //generate vector of landmarks within sensor_range of particle's location

  vector<LandmarkObs> close_landmarks;
  int num_landmarks = map_landmarks.landmark_list.size();
  for(int j = 0; j < num_landmarks; j++){
    Map::single_landmark_s lm = map_landmarks.landmark_list[j];
    double lm_dist = dist(pcle.x,pcle.y,lm.x_f,lm.y_f);
    if (lm_dist < sensor_range){
      LandmarkObs newlm;
      //newlm.id = lm.i_i; //should I use j instead?
      newlm.id = j; 
      newlm.x  = lm.x_f;
      newlm.y  = lm.y_f;
      close_landmarks.push_back(newlm);
    }
  }//for j over num_landmarks
  return close_landmarks;
}

vector<LandmarkObs> ParticleFilter::calc_trans_meas(Particle &pcle, const vector<LandmarkObs> &observations){
  //translate the observations from the vehicle's frame of reference to the map's coordinates

  vector<LandmarkObs> trans_measurements;
  for(int j = 0; j < num_meas; j++){
    LandmarkObs lobs = observations[j];
    LandmarkObs trans_measurement;
    trans_measurement.id = lobs.id;
    //does this need to change for non-standard coordinate orientation?  Doesn't seem like it
    trans_measurement.x = pcle.x + (cos(pcle.theta) * lobs.x) - (sin(pcle.theta) * lobs.y);
    trans_measurement.y = pcle.y + (sin(pcle.theta) * lobs.x) + (cos(pcle.theta) * lobs.y);
    trans_measurements.push_back(trans_measurement);
  }//for j over num_meas
  return trans_measurements;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  //cout<<"associating measurements with landmarks\n";
  for(int i=0; i<observations.size(); i++){
    //cout<<"  starting associating measurement "<<i<<"\n";
    double min_dist = 1e6;
    for(int j=0; j<predicted.size(); j++){
      //cout<<"    considering landmark "<<j<<"\n";
      double cur_dist = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
      if (cur_dist < min_dist){
        //cout<<"    updating min_dist from "<<min_dist<<" to "<<cur_dist<<"\n";
        observations[i].id = predicted[j].id;
        min_dist = cur_dist;
      }
    }//for j over predicted
    //cout<<"  finished associating measurement "<<i<<" with dist of "<<min_dist<<"\n";
  }//for i over observations
  //cout<<"  finished dataAssociation \n";
  return;
}

void ParticleFilter::calc_weight(Particle &pcle,std::vector<LandmarkObs> trans_measurements, const Map &map_landmarks, double std_landmark[]){
  //caluclate the weight for a particle, given its obesrvations and the associated landmarks from the map

  double pcle_prob = 1;
  //vector<int> associations;
  //vector<double> sense_x;
  //vector<double> sense_y;
  for(int j = 0; j < num_meas; j++){
    LandmarkObs& tm = trans_measurements[j];
    Map::single_landmark_s asc_lm = map_landmarks.landmark_list[tm.id];
    //associations.push_back(asc_lm.id_i);
    //sense_x.push_back(tm.x);
    //sense_y.push_back(tm.y);
    //particles[i] = SetAssociations(particles[i],associations,sense_x,sense_y);
    //cout<<"    calculating meas_prob with x's of "<<tm.x<<","<<asc_lm.x_f<<" and y's of "<<tm.y<<","<<asc_lm.y_f<<"\n";
    double meas_prob =
      //(1/(2*M_PI*std_landmark[0]*std_landmark[1])) * //ignore constant, since we will normalize anyway, and std_landmark is the same for all particles and measurements
      exp(
          -1 * (
                ((pow(tm.x - asc_lm.x_f,2))/(2*std_landmark[0]*std_landmark[0])) +
                ((pow(tm.y - asc_lm.y_f,2))/(2*std_landmark[1]*std_landmark[1]))
                )
          );
    //cout<<"    meas_prob is "<<meas_prob<<"\n";
    pcle_prob *= meas_prob;
  }
  pcle.weight = pcle_prob;
  //cout<<"  pcle_prob is "<<pcle_prob<<"\n";
  return;
}

int ParticleFilter::normalize_weights(){
  //normalize the weights, and return the index of the max weight, so that its associations can be set

  //cout<<"normalizing weights\n";
  double weight_sum = 0;
  double max_weight = 0;
  int max_index = -1;
  for(int i = 0; i< num_particles; i++){
    weight_sum += particles[i].weight;
    if (particles[i].weight > max_weight){
      max_weight = particles[i].weight;
      max_index = i;
    }
  }
  //cout<<"weight_sum is "<<weight_sum<<"\n";
  for(int i = 0; i< num_particles; i++){
    particles[i].weight /= weight_sum;
    weights[i] = particles[i].weight;
  }//for i over num_particles
  return max_index;
}

void ParticleFilter::prep_associations(int mi, double sensor_range,const Map &map_landmarks,const std::vector<LandmarkObs> &observations){
  //for particle with index mi, generate associations, sense_x, and sense_y vectors, and call set_associations

  //regenerate observation-landmark associations, since I didn't save them (and I think re-calcuating for a single particle is cheaper than saving for all particles)
  vector<LandmarkObs> close_landmarks = find_close_landmarks(sensor_range,particles[mi],map_landmarks);
  vector<LandmarkObs> trans_measurements = calc_trans_meas(particles[mi],observations);
  dataAssociation(close_landmarks,trans_measurements);

  //generate vectors for associations, sense_x, sense_y
  vector<int> associations;
  vector<double> sense_x;
  vector<double> sense_y;
  for(int j = 0; j < num_meas; j++){
    LandmarkObs& tm = trans_measurements[j];
    Map::single_landmark_s asc_lm = map_landmarks.landmark_list[tm.id];
    associations.push_back(asc_lm.id_i);
    sense_x.push_back(tm.x);
    sense_y.push_back(tm.y);
  }// for j over num_meas

  //call SetAssociations
  particles[mi] = SetAssociations(particles[mi],associations,sense_x,sense_y);
  return;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  //cout<<"starting updateWeights\n";
  num_meas = observations.size();
  for(int i = 0; i < num_particles; i++){
    Particle& pcle = particles[i];

    //find landmarks within sensor range
    //cout<<"finding landmarks within sensor range\n";
    vector<LandmarkObs> close_landmarks = find_close_landmarks(sensor_range,pcle,map_landmarks);

    //translate measurements to absolute coordinates for this particular particle
    //cout<<"translating measurements to absolute coordinates\n";
    vector<LandmarkObs> trans_measurements = calc_trans_meas(pcle,observations);

    //associate observations with closest landmark
    dataAssociation(close_landmarks,trans_measurements);

    //calculate weight for particle
    //cout<<"calculating particle's weight\n";
    calc_weight(pcle,trans_measurements,map_landmarks,std_landmark);

  }//for i over particles
  int mi = normalize_weights();
  prep_associations(mi,sensor_range,map_landmarks,observations);

  return;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  default_random_engine gen;
  discrete_distribution<int> disdist(weights.begin(), weights.end());
  vector<Particle> new_particles;
  for (int i = 0; i <num_particles; i++){
    new_particles.push_back(particles[disdist(gen)]);
  }
  particles = new_particles;
  return;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
