#ifndef __ICP_WRAPPER_H__
#define __ICP_WRAPPER_H__

#include <string>

double getInitialTransformationICP(const std::string& model_path, const std::string& template_path, const std::string& correspondences_path,
                                   Matrix &R, Matrix &t);

#endif //__ICP_WRAPPER_H__
