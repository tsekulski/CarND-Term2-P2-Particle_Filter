/*
 * particle_filter.cpp
 *
 *Implementation by T. Sekulski based on class structure and methods declaration by Tiffany Huang / MBRDNA / Udacity
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

#include "particle_filter.h"

using namespace std;

// Random engine - to be used in methods below
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 50;

	// Set the size of the particles & weights vectors
	particles.resize(num_particles);
	weights.resize(num_particles);

	// Standard deviations for x, y, and theta
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// Create normal distributions for x, y and theta
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; ++i) {

		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1.0;
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Standard deviations for x, y, and yaw
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	// Create normal distributions for x, y and theta - as white noise with mean = 0
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	for (int i = 0; i < num_particles; ++i) {

		//Predict position without Gaussian noise

		//Avoid division by zero
		if (fabs(yaw_rate) < 0.00001) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		else {
			double new_theta = particles[i].theta + yaw_rate*delta_t;
			particles[i].x += (velocity/yaw_rate)*(sin(new_theta) - sin(particles[i].theta));
			particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(new_theta));
			particles[i].theta = new_theta;
		}

		//Add Gaussian noise to predicted values
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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


	// ******** Iterate through all particles *******
	for (int i = 0; i < num_particles; ++i) {

		double x_part = particles[i].x;
		double y_part = particles[i].y;
		double theta = particles[i].theta;

		//Reset weight to 1.0 for further calculations
		particles[i].weight = 1.0;

		//Limit the list of landmarks to those within the sensor's range
		vector<LandmarkObs> landmarks_in_range;
		for (int k = 0; k < map_landmarks.landmark_list.size(); ++k){
			double temp_distance = dist(x_part, y_part, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
			if (temp_distance <= sensor_range) {
				landmarks_in_range.push_back(LandmarkObs{map_landmarks.landmark_list[k].id_i, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f});
			}
		}

		// ******** Iterate through all particle's observations ********
		for (int j = 0; j < observations.size(); j++){

			//1. Transform observations from the VEHICLE's coordinate system into the MAP's coordinate system
			double x_obs = observations[j].x;
			double y_obs = observations[j].y;

			// Transform to map x coordinate
			double x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);

			// Transform to map y coordinate
			double y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);

			//Append to transformed observations vector
			//particles[i].sense_x.push_back(x_map);
			//particles[i].sense_y.push_back(y_map);

			//2. Associate the transformed observations with the nearest landmark

			//Initialize distance & association with the first landmark
			double distance = dist(x_map, y_map, landmarks_in_range[0].x, landmarks_in_range[0].y);
			int association = 0;

			//Iterate through landmarks and check which one is nearest
			for (int k = 0; k < landmarks_in_range.size(); ++k){
				double temp_distance = dist(x_map, y_map, landmarks_in_range[k].x, landmarks_in_range[k].y);
				if (temp_distance <= distance) {
					distance = temp_distance;
					association = k;
				}
			}
			// Add association to the particle's associations list
			//particles[i].associations.push_back(association);

			//3. Calculate the partial weight related to the observation j

			// Define inputs
			double sig_x = std_landmark[0];
			double sig_y = std_landmark[1];
			x_obs = x_map;
			y_obs = y_map;
			double mu_x = landmarks_in_range[association].x;
			double mu_y = landmarks_in_range[association].y;

			// Calculate normalization term
			double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));

			// Calculate exponent
			double exponent = (pow(x_obs-mu_x, 2))/(2 * pow(sig_x, 2)) + (pow(y_obs-mu_y, 2))/(2 * pow(sig_y, 2));

			// Calculate weight using normalization terms and exponent
			double partial_weight = gauss_norm * exp(-exponent);

			//4. To get the final weight just multiply all the calculated measurement probabilities together
			particles[i].weight *= partial_weight;

		} // ******** Finished iterating through all particle's observations ********

		//Update weights vector
		weights[i] = particles[i].weight;

		//Clear particles data on observations and associations to prepare for the next update run
		//particles[i].sense_x.clear();
		//particles[i].sense_y.clear();
		//particles[i].associations.clear();

	} // ******** Finished iterating through all particles *******

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	/* Python implementation - resampling wheel
	p3 = []
	index = int(random.random()*N)
	beta = 0.0
	mw = max(w)
	for i in range(N):
    	beta += random.random() * 2.0 * mw
    	while beta > w[index]:
        	beta -= w[index]
        	index = (index + 1) % N
    	p3.append(p[index])
	p = p3
	*/

	//Vector to hold resampled particles
	vector<Particle> resampled_particles(num_particles);

	//Random particle distribution with probability according to their weights
	discrete_distribution<int> disc_dist(weights.begin(), weights.end());

	//Draw random particles according to probability distribution set above
	for (int i = 0; i < num_particles; ++i) {
		resampled_particles[i] = particles[disc_dist(gen)];
	}

	//Update the particles vector
	particles = resampled_particles;
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
