#ifndef ABSTRACTZCLUSTER_H
#define ABSTRACTZCLUSTER_H

#include <string>
#include <ginac/ginac.h>

#include "MyLambdaPolynomial.h"

/// abstract class for partition function on a cluster
class AbstractZCluster {
protected:
    std::string Type;

    virtual void EvaluateZ() = 0; /// evaluate Z on the cluster

public:
    /// constructor
    AbstractZCluster(const std::string& type) : Type(type) {}

    /// destructor
    virtual ~AbstractZCluster() = default;

    virtual MyLambdaPolynomial<GiNaC::numeric> ComputeLambdaPolynomial() = 0;

    virtual void PrintZ() const = 0;

    std::string GetType() const { return this->Type; }
};

#endif // ABSTRACTZCLUSTER_H
