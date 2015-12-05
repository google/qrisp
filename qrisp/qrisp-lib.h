#include "dataset-utils.h"
#include "qrisp/proto/config.pb.h"
#include "rna-structure.h"
#include "utils.h"

using namespace std;

namespace qrisp {

void StartPrediction(const Dataset& dataset, const QRSPModel& model,
                     const Configuration& config);

void StartTraining(const Dataset& training_set, const Dataset& holdout_set,
                   const QRSPModel& model, const Configuration& config);

}  // namespace qrisp
