#!/bin/bash

code=`git log --pretty=format:'%h' -n 1`
echo "#ifndef GIT_VERSION_H
#define GIT_VERSION_H
#include <string>
static std::string git_version = \"$code\";
#endif" > source/git_version.h

