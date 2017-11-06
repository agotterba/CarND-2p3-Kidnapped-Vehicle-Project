/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};



class ParticleFilter {
	
  // Number of particles to draw
  int num_particles; 

  // Number of measurements for current timestep
  int num_meas;
	
  // Flag, if filter is initialized
  bool is_initialized;
	
  // Vector of weights of all particles
  std::vector<double> weights;

  //limit for yaw_rate to avoid dividing by zero
  double yaw_rate_limit;

  //random generators for prediction
  //default_random_engine pred_gen;
  //normal_distribution<double> pred_dist_x(0, std[0]);
  //normal_distribution<double> pred_dist_y(0, std[1]);
  //normal_distribution<double> pred_dist_theta(0,std[2]);

public:
	
  // Set of current particles
  std::vector<Particle> particles;

  // Constructor
  // @param num_particles Number of particles
  ParticleFilter() : num_particles(0), is_initialized(false) {}

  // Destructor
  ~ParticleFilter() {}

  /**
   * init Initializes particle filter by initializing particles to Gaussian
   *   distribution around first position and all the weights to 1.
   * @param x Initial x position [m] (simulated estimate from GPS)
   * @param y Initial y position [m]
   * @param theta Initial orientation [rad]
   * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
   *   standard deviation of yaw [rad]]
   */
  void init(double x, double y, double theta, double std[]);
  
  /**
   * prediction Predicts the state for the next time step
   *   using the process model.
   * @param delta_t Time between time step t and t+1 in measurements [s]
   * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
   *   standard deviation of yaw [rad]]
   * @param velocity Velocity of car from t to t+1 [m/s]
   * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
   */
  void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
  
  /**
   * dataAssociation Finds which observations correspond to which landmarks (likely by using
   *   a nearest-neighbors data association).
   * @param predicted Vector of predicted landmark observations
   * @param observations Vector of landmark observations
   */
  void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
  
  /**
   * updateWeights Updates the weights for each particle based on the likelihood of the 
   *   observed measurements. 
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
   * @param observations Vector of landmark observations
   * @param map Map class containing map landmarks
   */
  void updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations,
                     const Map &map_landmarks);
  
  /**
   * resample Resamples from the updated set of particles to form
   *   the new set of particles.
   */
  void resample();
  
  /*
   * Set a particles list of associations, along with the associations calculated world x,y coordinates
   * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
   */
  Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);
  
  std::string getAssociations(Particle best);
  std::string getSenseX(Particle best);
  std::string getSenseY(Particle best);

  //generate vector of landmarks within sensor_range of particle's location
  std::vector<LandmarkObs> find_close_landmarks(double sensor_range, Particle &pcle, const Map &map_landmarks);

  //translate the observations from the vehicle's frame of reference to the map's coordinates
  std::vector<LandmarkObs> calc_trans_meas(Particle &pcle, const std::vector<LandmarkObs> &observations);

  //caluclate the weight for a particle, given its obesrvations and the associated landmarks from the map
  void calc_weight(Particle &pcle,std::vector<LandmarkObs> trans_measurements, const Map &map_landmarks, double std_landmark[]);

  //normalize the weights, and return the index of the max weight, so that its associations can be set
  int normalize_weights();

  //for particle with index mi, generate associations, sense_x, and sense_y vectors, and call set_associations
  void prep_associations(int mi,double sensor_range,const Map &map_landmarks,const std::vector<LandmarkObs> &observations);
  
  /**
   * initialized Returns whether particle filter is initialized yet or not.
   */
  const bool initialized() const {
    return is_initialized;
  }
};

#endif /* PARTICLE_FILTER_H_ */
