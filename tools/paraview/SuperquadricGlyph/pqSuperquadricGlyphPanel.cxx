

#include "pqSuperquadricGlyphPanel.h"


pqSuperquadricGlyphPanel::pqSuperquadricGlyphPanel(pqProxy *proxy, QWidget *p)
  : Superclass(proxy, p)
{
  // -------------------------------------------------------------------------
  // Automatically enable roundness widgets if roundness should not be set individually per object.

  this->FixedThetaPhiRoundnessWidget = this->findChild<QCheckBox*>("FixedThetaPhiRoundness");
  if (! this->FixedThetaPhiRoundnessWidget)
  {
    qWarning() << "Failed to locate FixedThetaPhiRoundness widget.";
    return;
  }

  this->ThetaRoundnessWidget = this->findChild<QWidget*>("ThetaRoundness");
  if (this->ThetaRoundnessWidget)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->ThetaRoundnessWidget, SLOT(setEnabled(bool)));

  this->ThetaRoundnessLabel = this->findChild<QLabel*>("_labelForThetaRoundness");
  if (this->ThetaRoundnessLabel)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->ThetaRoundnessLabel, SLOT(setEnabled(bool)));


  this->PhiRoundnessWidget = this->findChild<QWidget*>("PhiRoundness");
  if (this->PhiRoundnessWidget)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->PhiRoundnessWidget, SLOT(setEnabled(bool)));

  this->PhiRoundnessLabel = this->findChild<QLabel*>("_labelForPhiRoundness");
  if (this->PhiRoundnessLabel)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->PhiRoundnessLabel, SLOT(setEnabled(bool)));

  this->FixedThetaPhiRoundnessWidget->toggle();
  this->FixedThetaPhiRoundnessWidget->toggle();

  // -------------------------------------------------------------------------

  this->HalfAxisVectorWidget = this->findChild<QCheckBox*>("HalfAxisVector");
  if (! this->HalfAxisVectorWidget)
  {
    qWarning() << "Failed to locate HalfAxisVector widget.";
    return;
  }

  this->HalfAxisVectorArrayWidget = this->findChild<QWidget*>("HalfAxisVectorArray");
  if (this->HalfAxisVectorArrayWidget)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisVectorArrayWidget, SLOT(setEnabled(bool)));

  this->HalfAxisVectorArrayLabel = this->findChild<QLabel*>("_labelForHalfAxisVectorArray");
  if (this->HalfAxisVectorArrayLabel)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisVectorArrayLabel, SLOT(setEnabled(bool)));


  this->HalfAxisXArrayWidget = this->findChild<QWidget*>("HalfAxisXArray");
  if (this->HalfAxisXArrayWidget)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisXArrayWidget, SLOT(setDisabled(bool)));

  this->HalfAxisXArrayLabel = this->findChild<QLabel*>("_labelForHalfAxisXArray");
  if (this->HalfAxisXArrayLabel)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisXArrayLabel, SLOT(setDisabled(bool)));

  this->HalfAxisYArrayWidget = this->findChild<QWidget*>("HalfAxisYArray");
  if (this->HalfAxisYArrayWidget)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisYArrayWidget, SLOT(setDisabled(bool)));

  this->HalfAxisYArrayLabel = this->findChild<QLabel*>("_labelForHalfAxisYArray");
  if (this->HalfAxisYArrayLabel)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisYArrayLabel, SLOT(setDisabled(bool)));

  this->HalfAxisZArrayWidget = this->findChild<QWidget*>("HalfAxisZArray");
  if (this->HalfAxisZArrayWidget)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisZArrayWidget, SLOT(setDisabled(bool)));

  this->HalfAxisZArrayLabel = this->findChild<QLabel*>("_labelForHalfAxisZArray");
  if (this->HalfAxisZArrayLabel)
    QObject::connect(this->HalfAxisVectorWidget, SIGNAL(toggled(bool)),
                     this->HalfAxisZArrayLabel, SLOT(setDisabled(bool)));

  this->HalfAxisVectorWidget->toggle();
  this->HalfAxisVectorWidget->toggle();

  // -------------------------------------------------------------------------

  this->QuatOrientationWidget = this->findChild<QCheckBox*>("QuatOrientation");
  if (! this->QuatOrientationWidget)
  {
    qWarning() << "Failed to locate QuatOrientation widget.";
    return;
  }

  this->Quat1ArrayWidget = this->findChild<QWidget*>("Quat1Array");
  if (this->Quat1ArrayWidget)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat1ArrayWidget, SLOT(setEnabled(bool)));

  this->Quat1ArrayLabel = this->findChild<QLabel*>("_labelForQuat1Array");
  if (this->Quat1ArrayLabel)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat1ArrayLabel, SLOT(setEnabled(bool)));

  this->Quat2ArrayWidget = this->findChild<QWidget*>("Quat2Array");
  if (this->Quat2ArrayWidget)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat2ArrayWidget, SLOT(setEnabled(bool)));

  this->Quat2ArrayLabel = this->findChild<QLabel*>("_labelForQuat2Array");
  if (this->Quat2ArrayLabel)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat2ArrayLabel, SLOT(setEnabled(bool)));

  this->Quat3ArrayWidget = this->findChild<QWidget*>("Quat3Array");
  if (this->Quat3ArrayWidget)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat3ArrayWidget, SLOT(setEnabled(bool)));

  this->Quat3ArrayLabel = this->findChild<QLabel*>("_labelForQuat3Array");
  if (this->Quat3ArrayLabel)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat3ArrayLabel, SLOT(setEnabled(bool)));

  this->Quat4ArrayWidget = this->findChild<QWidget*>("Quat4Array");
  if (this->Quat4ArrayWidget)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat4ArrayWidget, SLOT(setEnabled(bool)));

  this->Quat4ArrayLabel = this->findChild<QLabel*>("_labelForQuat4Array");
  if (this->Quat4ArrayLabel)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->Quat4ArrayLabel, SLOT(setEnabled(bool)));

  this->RotationMatrixArrayWidget = this->findChild<QWidget*>("RotationMatrixArray");
  if (this->RotationMatrixArrayWidget)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->RotationMatrixArrayWidget, SLOT(setDisabled(bool)));

  this->RotationMatrixArrayLabel = this->findChild<QLabel*>("_labelForRotationMatrixArray");
  if (this->RotationMatrixArrayLabel)
    QObject::connect(this->QuatOrientationWidget, SIGNAL(toggled(bool)),
                     this->RotationMatrixArrayLabel, SLOT(setDisabled(bool)));

  this->QuatOrientationWidget->toggle();
  this->QuatOrientationWidget->toggle();

  // -------------------------------------------------------------------------


  return;
}
