/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2023, Michal Minařík
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
 *********************************************************************/

#ifndef OMPL_GEOMETRIC_PLANNERS_RRT_RRT_LIB_
#define OMPL_GEOMETRIC_PLANNERS_RRT_RRT_LIB_

#include "ompl/datastructures/NearestNeighbors.h"
#include "ompl/geometric/planners/PlannerIncludes.h"

namespace ompl
{
    namespace geometric
    {
        /**
           @anchor gRRTLIB
        */

        /** \brief RRT-LIB */
        class RRTLIB : public base::Planner
        {
        public:
            /** \brief Constructor */
            RRTLIB(const base::SpaceInformationPtr &si, bool addIntermediateStates = false);

            ~RRTLIB() override;

            void getPlannerData(base::PlannerData &data) const override;

            base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc) override;

            void clear() override;

            /** \brief Set the goal bias

                In the process of randomly selecting states in
                the state space to attempt to go towards, the
                algorithm may in fact choose the actual goal state, if
                it knows it, with some probability. This probability
                is a real number between 0.0 and 1.0; its value should
                usually be around 0.05 and should not be too large. It
                is probably a good idea to use the default value. */
            void setGoalBias(double goalBias)
            {
                goalBias_ = goalBias;
            }

            /** \brief Get the goal bias the planner is using */
            double getGoalBias() const
            {
                return goalBias_;
            }

            /** \brief Set the path bias */
            void setPathBias(double pathBias)
            {
                pathBias_ = pathBias;
            }

            /** \brief Get the path bias the planner is using */
            double getPathBias() const
            {
                return pathBias_;
            }

            /** \brief Set the safe distance around start and goal states
                       where inhibited regions are not created bias */
            void setSafeDistance(double safeDistance)
            {
                safeDistance_ = safeDistance;
            }

            /** \brief Get the safe distance the planner is using */
            double getSafeDistance() const
            {
                return safeDistance_;
            }

            /** \brief Set the inhibited region redius */
            void setInhibitedRadius(double inhibitedRadius)
            {
                inhibitedRadius_ = inhibitedRadius;
            }

            /** \brief Get the inhibited region redius the planner is using */
            double getInhibitedRadius() const
            {
                return inhibitedRadius_;
            }

            /** \brief Set the guiding redius */
            void setGuidingRadius(double guidingRadius)
            {
                guidingRadius_ = guidingRadius;
            }

            /** \brief Get the guiding redius the planner is using */
            double getGuidingRadius() const
            {
                return guidingRadius_;
            }

            /** \brief Set the diversity patience */
            void setDiversityPatience(int diversityPatience)
            {
                diversityPatience_ = diversityPatience;
            }

            /** \brief Get the diversity patience the planner is using */
            int getDiversityPatience() const
            {
                return diversityPatience_;
            }
            /** \brief Set the guiding paths directory */
            void setGuidingPathsDirectory(std::string guidingPathsDirectory)
            {
                guidingPathsDirectory_ = guidingPathsDirectory;
            }

            /** \brief Get the guiding paths directory the planner is using */
            std::string getGuidingPathsDirectory() const
            {
                return guidingPathsDirectory_;
            }

            /** \brief Set the guiding paths indices */
            void setGuidingPathsIndices(std::string guidingPathsIndices)
            {
                guidingPathsIndices_ = guidingPathsIndices;
            }

            /** \brief Get the guiding paths indices the planner is using */
            std::string getGuidingPathsIndices() const
            {
                return guidingPathsIndices_;
            }

            void setGuidingPathsTransformation(const double T[4][4])
            {
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        T_[i][j] = T[i][j];
                    }
                }
            }

            /** \brief Return true if the intermediate states generated along
                       motions are to be added to the tree itself */
            bool getIntermediateStates() const
            {
                return addIntermediateStates_;
            }

            /* \brief Specify whether the intermediate states generated along
                      motions are to be added to the tree itself */
            void setIntermediateStates(bool addIntermediateStates)
            {
                addIntermediateStates_ = addIntermediateStates;
            }

            /* \brief Return true if the paths will be saved */
            bool getSaveResults() const
            {
                return saveResults_;
            }

            /** \brief Specify whether the paths should be saved */
            void setSaveResults(bool saveResults)
            {
                saveResults_ = saveResults;
            }

            /* \brief Return the output folder */
            std::string getOutput() const
            {
                return output_;
            }

            /** \brief Specify output folder */
            void setOutput(std::string output)
            {
                output_ = output;
            }

            /** \brief Return used mode */
            int getMode() const
            {
                return mode_;
            }

            /** \brief Specify run mode
                        0 - generate guiding paths and plan
                        1 - load guiding paths and plan
                        2 - generate guiding paths only */
            void setMode(int mode)
            {
                mode_ = mode;
            }

            /** \brief Set the range the planner is supposed to use.

                This parameter greatly influences the runtime of the
                algorithm. It represents the maximum length of a
                motion to be added in the tree of motions. */
            void setRange(double distance)
            {
                maxDistance_ = distance;
            }

            /** \brief Get the range the planner is using */
            double getRange() const
            {
                return maxDistance_;
            }

            /** \brief Set the names of files used

                The filenames are specified only to be saved in the
                info file. Planner behavior is not modified. */
            void setFilenames(const std::string &map, const std::string &guidingObject, const std::string &object, const std::string &correspondences)
            {
                fileMap_ = map;
                fileGuidingObject_ = guidingObject;
                fileObject_ = object;
                fileCorrespondences_ = correspondences;
            }

            /** \brief Set the problem setup

                The parameters are specified only to be saved in the
                info file. Planner behavior is not modified. */
            void setSetup(const std::string &startConfiguration, const std::string &goalConfiguration, const std::string &spaceBounds)
            {
                startConfigurationStr_ = startConfiguration;
                goalConfigurationStr_ = goalConfiguration;
                spaceBoundsStr_ = spaceBounds;
            }

            /** \brief Set a different nearest neighbors datastructure */
            template <template <typename T> class NN>
            void setNearestNeighbors()
            {
                if (nn_ && nn_->size() != 0)
                    OMPL_WARN("Calling setNearestNeighbors will clear all states.");
                clear();
                nn_ = std::make_shared<NN<Motion *>>();
                setup();
            }

            void setup() override;

        protected:
            /** \brief Representation of a motion

                This only contains pointers to parent motions as we
                only need to go backwards in the tree. */
            class Motion
            {
            public:
                Motion() = default;

                /** \brief Constructor that allocates memory for the state */
                Motion(const base::SpaceInformationPtr &si) : state(si->allocState())
                {
                }

                ~Motion() = default;

                /** \brief The state contained by the motion */
                base::State *state{nullptr};

                /** \brief The parent motion in the exploration tree */
                Motion *parent{nullptr};
            };

            /** \brief Representation of a inhibited region */
            class InhibitedRegion
            {
            public:
                InhibitedRegion() = default;

                /** \brief Constructor that allocates memory for the state */
                InhibitedRegion(const base::SpaceInformationPtr &si, const int &pathIndex, const int &regionIndex) : state(si->allocState()), pathIndex(pathIndex), regionIndex(regionIndex)
                {
                }

                ~InhibitedRegion() = default;

                /** \brief Region center */
                base::State *state{nullptr};

                /** \brief Index of guiding path this region is a part of */
                size_t pathIndex;

                /** \brief Index of this region in the parent guiding path */
                size_t regionIndex;
            };

            /** \brief Try to find a path along specified guiding paths

                Inhibited regions are not pass as a parameter, member variable
                inhibitedNN_ is used instead. This is to enable more efficient
                control of the memory used by the kd-tree  llowing for NN
                search within inhibited regions. */
            base::PlannerStatus findPath(const base::PlannerTerminationCondition &ptc, const std::vector<ompl::geometric::PathGeometricPtr> &guidingPaths);

            /** \brief Saves information about the run and paths passed in `paths` to a specified folder */
            void saveResults(const std::vector<ompl::geometric::PathGeometricPtr> &paths) const;

            /** \brief Loads a guiding path at `pathFilename` and adds it to `guidingPaths` */
            bool loadGuidingPath(const std::string &pathFilename, std::vector<ompl::geometric::PathGeometricPtr> &guidingPaths);

            /** \brief Calculates the probability of entering an inhibited region  */
            inline double calculateP_(const size_t &i, const size_t &j, const std::vector<std::vector<size_t>> &attempts, const size_t &attempts_sum) const;

            /** \brief Calculates the similarity of two paths */
            double getSimilarity(const std::shared_ptr<ompl::geometric::PathGeometric> &a, const std::shared_ptr<ompl::geometric::PathGeometric> &b) const;

            /** \brief Calculates the distinctness between `new_path` and member variable `guidingPaths_` */
            double pathDistinctness(const std::shared_ptr<PathGeometric> &new_path) const;

            /** \brief Free the memory allocated by this planner */
            void freeMemory();
            void freeNN_();
            void freeInhibitedNN_();

            /** \brief Compute distance between motions (actually distance between contained states) */
            double distanceFunction(const Motion *a, const Motion *b) const
            {
                return si_->distance(a->state, b->state);
            }

            /** \brief Compute distance between motions (actually distance between contained states) */
            double distanceFunction(const InhibitedRegion *a, const InhibitedRegion *b) const
            {
                return si_->distance(a->state, b->state);
            }

            /** \brief State sampler */
            base::StateSamplerPtr sampler_;

            /** \brief A nearest-neighbors datastructure containing the tree of motions */
            std::shared_ptr<NearestNeighbors<Motion *>> nn_;

            std::shared_ptr<NearestNeighbors<InhibitedRegion *>> inhibitedNN_;

            /** \brief 0 - generate guiding paths prior to the planning
                       1 - load guiding paths prior to the planning
                       2 - generate guiding paths only */
            int mode_{0};

            /** \brief The fraction of time the goal is picked as the state to expand towards (if such a state is available) */
            double goalBias_{.05};

            /** \brief The fraction of time a state along the guiding paths is picked as the state to expand towards (if such a state is available) */
            double pathBias_{.80};

            /** \brief The radius around the start and the goal configuration in which configurations are not added to the inhibited regions */
            double safeDistance_{2.00};

            /** \brief The radius of inhibited regions */
            double inhibitedRadius_{3.00};

            /** \brief The radius around the guiding path configurations where new configurations are sampled */
            double guidingRadius_{0.50};

            /** \brief Number of non-distinct paths needed to find before the search is stopped */
            int diversityPatience_{10};

            /** \brief Directory containing the guiding paths to be used */
            std::string guidingPathsDirectory_{""};

            /** \brief Indices of the guiding paths to be used */
            std::string guidingPathsIndices_{""};

            /** \brief Object filenames (for saving purposes only) */
            std::string fileMap_{""}, fileGuidingObject_{""}, fileObject_{""}, fileCorrespondences_{""};

            /** \brief Problem setup (for saving purposes only) */
            std::string startConfigurationStr_{""}, goalConfigurationStr_{""}, spaceBoundsStr_{""};

            /** \brief The maximum length of a motion to be added to a tree */
            double maxDistance_{0.};

            /** \brief Flag indicating whether intermediate states are added to the built tree of motions */
            bool addIntermediateStates_;

            /** \brief Flag indicating whether paths will be saved */
            bool saveResults_{false};

            std::string output_{""};

            /** \brief The random number generator */
            RNG rng_;

            /** \brief The most recent goal motion.  Used for PlannerData computation */
            Motion *lastGoalMotion_{nullptr};

            /** \brief The most recent found path.  Used for inhibited region building */
            std::shared_ptr<PathGeometric> lastFoundPath_{nullptr};

            /** \brief Guiding paths found during the initial search */
            std::vector<ompl::geometric::PathGeometricPtr> guidingPaths_{};

            double T_[4][4]{{1, 0, 0, 0},
                            {0, 1, 0, 0},
                            {0, 0, 1, 0},
                            {0, 0, 0, 1}};
        };
    }
}

#endif
