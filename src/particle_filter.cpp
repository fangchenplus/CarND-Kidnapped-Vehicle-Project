/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;


double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  // TODO: Set the number of particles
  
  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  // Set standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2]; 

  // Creates normal (Gaussian) distributions for x, y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  for (int i = 0; i < num_particles; ++i) {
    double sample_x, sample_y, sample_theta;
    
    // Sample from these normal distributions like this: 
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);   
     
    // Print your samples to the terminal.
    // std::cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " 
              // << sample_theta << std::endl;
			  
    Particle p;
    p.id = i;
    p.x = sample_x;
    p.y = sample_y;
    p.theta = sample_theta;
    p.weight = 1.0;
    
    particles.push_back(p);
    weights.push_back(p.weight);
  }
  
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  // Set standard deviations for x, y, and theta
  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2]; 

  // Creates normal (Gaussian) distributions for x, y and theta
  normal_distribution<double> dist_x(0.0, std_x);
  normal_distribution<double> dist_y(0.0, std_y);
  normal_distribution<double> dist_theta(0.0, std_theta);
  
  // std::cout << "prediction " << yaw_rate << std::endl;   

  for (int i = 0; i < num_particles; ++i) {
    double sample_x, sample_y, sample_theta;
    
    // Sample from these normal distributions like this: 
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);   
    
    double theta0, theta_new;
    theta0 = particles[i].theta;
    theta_new = theta0 + yaw_rate*delta_t;
    
    if (yaw_rate == 0){
      particles[i].x += velocity*delta_t * cos(theta0) ;
      particles[i].y += velocity*delta_t * sin(theta0) ;
    }
    else {
      particles[i].x += velocity/yaw_rate* ( sin(theta_new) - sin(theta0) ) ;
      particles[i].y += velocity/yaw_rate* ( cos(theta0) - cos(theta_new) ) ;
    }
  
    particles[i].theta = theta_new ;
    
    particles[i].x += sample_x;
    particles[i].y += sample_y;
    particles[i].theta += sample_theta;      

  }   

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // int id_tmp; 
  double dist_min;
  double dist_tmp;
  for (unsigned int i = 0; i < observations.size(); i++) {
    dist_min = 1000.0;
    for (unsigned int j = 0; j < predicted.size(); j++) {
      dist_tmp = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
      if (dist_tmp < dist_min) {
        observations[i].id = predicted[j].id;
        dist_min = dist_tmp;
      }
    }
  // std::cout << "dataAssociation " << observations[i].id << std::endl;    
  }
    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for (int i = 0; i < num_particles; ++i) {
    double x_part, y_part, theta_part;
    x_part = particles[i].x;
    y_part = particles[i].y;
    theta_part = particles[i].theta;
 
    double x_map, y_map;   

    LandmarkObs single_landmark_tmp;
    
    // calculate the predicted observations
    vector<LandmarkObs> predicted;    
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      x_map = map_landmarks.landmark_list[j].x_f;      
      y_map = map_landmarks.landmark_list[j].y_f; 
      if (dist(x_part,y_part,x_map,y_map) <= sensor_range) {
        single_landmark_tmp.id = map_landmarks.landmark_list[j].id_i;
        single_landmark_tmp.x = x_map;
        single_landmark_tmp.y = y_map;
        predicted.push_back(single_landmark_tmp);
      }
    }
    
    
    // calculate the real observations in map coordinates
    vector<LandmarkObs> observations_map;
    for (unsigned int j = 0; j < observations.size(); j++) {
      
      // int id_obs = observations[j].id;
      double x_obs = observations[j].x;
      double y_obs = observations[j].y;
      
      // transform to map coordinate
      x_map = x_part + (cos(theta_part) * x_obs) - (sin(theta_part) * y_obs);
      y_map = y_part + (sin(theta_part) * x_obs) + (cos(theta_part) * y_obs);    
      
      single_landmark_tmp.id = -1;
      single_landmark_tmp.x = x_map;
      single_landmark_tmp.y = y_map;      
      observations_map.push_back(single_landmark_tmp);   
      
    }
  
    dataAssociation(predicted, observations_map); 
    
    double weight_tmp = 1.0;
    bool found = false;
    // double id_tmp;
    for (unsigned int j = 0; j < observations_map.size(); j++) {
      found = false;
      for (unsigned int k = 0; k < predicted.size(); k++) { 
        if ( predicted[k].id == observations_map[j].id){
          weight_tmp *= multiv_prob(std_landmark[0], std_landmark[1], observations_map[j].x, observations_map[j].y, predicted[k].x, predicted[k].y);
          found = true;
        } 
      // weight_tmp *= 0.0000001; // cannot find corresonding observation  
      }
      if (found == false) {
        weight_tmp *= 0.0000001; // cannot find corresonding observation 
      }

    }    
    // std::cout << "updateWeights " << weight_tmp << " " << predicted.size() << " " << observations_map.size() << std::endl;     
    particles[i].weight = weight_tmp;  
    weights[i] = weight_tmp;
  }     

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;

  std::discrete_distribution<int> dist_p(weights.begin(),weights.end());

  std::vector<Particle> new_particles;
  int sample_p;
  for (int i = 0; i < num_particles; ++i) {
    sample_p = dist_p(gen);    
    new_particles.push_back(particles[sample_p]);
  }
  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}