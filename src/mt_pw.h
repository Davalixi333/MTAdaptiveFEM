#ifndef _EM_PW_H_
#define _EM_PW_H_ 1

#include "mt_defs.h"

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include <memory>

class PlaneWaveFunction {
public:
  void initialize(const std::vector<Layer> &layered_model, double freq) {
    int i;
    Complex gmogp, rjexp, Atop;
    std::vector<Complex> kk, Rp, expmgh;

    nlay_ = (int)layered_model.size();

    ztop_.resize(nlay_);
    lsig_.resize(nlay_);

    for (i = 0; i < nlay_; ++i) {
      ztop_[i] = layered_model[i].first;
      lsig_[i] = layered_model[i].second;
    }

    omega_ = 2 * dealii::numbers::PI * freq;

    kk.resize(nlay_ + 1);
    Rp.resize(nlay_ + 1);
    expmgh.resize(nlay_ + 1);

    for (i = 1; i <= nlay_; ++i) {
      kk[i] = std::sqrt(-II * omega_ * MU * lsig_[i - 1]);
    }
    kk[0] = std::sqrt(-II * omega_ * MU * lsig_[0]);

    std::fill(expmgh.begin(), expmgh.end(), Complex(0.0));
    for (i = 1; i <= nlay_ - 1; ++i) {
      expmgh[i] = std::exp(-kk[i] * (ztop_[i] - ztop_[i - 1]));
    }

    std::fill(Rp.begin(), Rp.end(), Complex(0.0));
    for (i = nlay_ - 1; i >= 1; --i) {
      gmogp = (-II * omega_ * MU * (lsig_[i - 1] - lsig_[i])) / std::pow(kk[i] + kk[i + 1], 2);
      rjexp = Rp[i + 1] * expmgh[i + 1];
      Rp[i] = (gmogp + rjexp) / (1.0 + gmogp * rjexp) * expmgh[i];
    }

    a.resize(nlay_ + 1);
    b.resize(nlay_ + 1);
    c.resize(nlay_ + 1);
    d.resize(nlay_ + 1);
    std::fill(a.begin(), a.end(), Complex(0.0));
    std::fill(b.begin(), b.end(), Complex(0.0));
    std::fill(c.begin(), c.end(), Complex(0.0));
    std::fill(d.begin(), d.end(), Complex(0.0));
    a[0] = Complex(1.0);
    c[0] = Complex(1.0);

    for (i = 1; i <= nlay_; ++i) {
      Atop = a[i - 1] + b[i - 1] * expmgh[i - 1];
      b[i] = Atop / (1.0 + Rp[i] * expmgh[i]);
      a[i] = b[i] * Rp[i];

      Atop = c[i - 1] + d[i - 1] * expmgh[i - 1];
      d[i] = Atop / (1.0 - Rp[i] * expmgh[i]);
      c[i] = -d[i] * Rp[i];
    }
  }

  void calculate_ehfield(double z, Complex eh[4]) const {
    int i, ilay;
    Complex kk, expp, expm, exppdz, expmdz;

    ilay = 0;
    for (i = nlay_ - 1; i >= 0; --i) {
      if (z >= ztop_[i]) {
        ilay = i;
        break;
      }
    }

    kk = std::sqrt(-II * omega_ * MU * lsig_[ilay]);

    expp = Complex(0.0);
    expm = Complex(0.0);
    exppdz = Complex(0.0);
    expmdz = Complex(0.0);

    expm = std::exp(-kk * (z - ztop_[ilay]));
    expmdz = -kk * expm;
    if (ilay != (nlay_ - 1)) {
      expp = std::exp(kk * (z - ztop_[ilay + 1]));
      exppdz = kk * expp;
    }

    eh[0] = a[ilay + 1] * expp + b[ilay + 1] * expm;
    eh[1] = (a[ilay + 1] * exppdz + b[ilay + 1] * expmdz) / (II * MU * omega_);
    eh[2] = (c[ilay + 1] * exppdz + d[ilay + 1] * expmdz) / lsig_[ilay];
    eh[3] = c[ilay + 1] * expp + d[ilay + 1] * expm;

    if (std::abs(eh[0]) < 1E-80) {
      eh[0] = 0.0;
    }
    if (std::abs(eh[1]) < 1E-80) {
      eh[1] = 0.0;
    }
    if (std::abs(eh[2]) < 1E-80) {
      eh[2] = 0.0;
    }
    if (std::abs(eh[3]) < 1E-80) {
      eh[3] = 0.0;
    }
  }

private:
  int nlay_;
  double omega_;
  std::vector<double> ztop_, lsig_;
  std::vector<Complex> a, b, c, d;
};

class BoundaryFunction : public dealii::Function<2, Complex> {
public:
  BoundaryFunction() : dealii::Function<2, Complex>(1) {}

  void reinit(const BoundingBox &bbox, const std::vector<Layer> *layered_models, double freq, int mode) {
    int c;

    this->bbox = bbox;
    this->mode = mode;

    for (c = 0; c < 2; ++c) {
      pwfs[c].reset(new PlaneWaveFunction);
      pwfs[c]->initialize(layered_models[c], freq);
    }
  }

  virtual Complex value(const Point &p, unsigned int component) const override {
    int c;
    Complex eh[2][4], v;
    double length, coeffs[2];

    for (c = 0; c < 2; ++c) {
      pwfs[c]->calculate_ehfield(p[1], eh[c]);
    }

    length = (bbox.second[0] - bbox.first[0]);

    coeffs[Left] = (bbox.second[0] - p[0])  / length;
    coeffs[Right] = (p[0] - bbox.first[0]) / length;

    v = 0.0;
    if (mode == TE) {
      for (c = 0; c < 2; ++c) {
        v += coeffs[c] * eh[c][0];
      }
    } else {
      for (c = 0; c < 2; ++c) {
        v += coeffs[c] * eh[c][3];
      }
    }

    (void)component;

    return std::conj(v);
  }

private:
  int mode;
  BoundingBox bbox;
  std::shared_ptr<PlaneWaveFunction> pwfs[2];
};

#endif
