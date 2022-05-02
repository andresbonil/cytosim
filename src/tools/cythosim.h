#ifndef REPORT_PYTHON_H
#define REPORT_PYTHON_H

//#include "simul.h"
#include "sim_thread.h"
#include "messages.h"
#include "glossary.h"
#include "exceptions.h"
#include "print_color.h"
#include "filepath.h"
#include "splash.h"
#include "tictoc.h"
#include <csignal>
#include "unistd.h"
#include "python_frame.h"
#include <functional>
namespace py = pybind11;

void bar(void);
class SimThread;
#endif
