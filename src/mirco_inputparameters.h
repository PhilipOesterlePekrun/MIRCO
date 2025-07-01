#ifndef SRC_INPUTPARAMETERS_H_
#define SRC_INPUTPARAMETERS_H_

#include <Teuchos_ParameterList.hpp>  //# these are for the second ctor. though if we should move away from teuchos, we can change the SetParameters signature
#include <Teuchos_RCP.hpp>            //#
#include <Teuchos_XMLParameterListHelpers.hpp>  //#
#include <string>

namespace MIRCO
{
  /**
   * @brief This class stores the input parameters and can assign them from the .xml or RMG
   *
   */
  class InputParameters
  {
   public:
    InputParameters(const std::string& inputFileName)
    {
      Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList());
      Teuchos::updateParametersFromXmlFile(inputFileName, paramList.ptr());
      SetParameters(paramList);
    }
    InputParameters(Teuchos::RCP<Teuchos::ParameterList> paramList) { SetParameters(paramList); }



   private:
    void SetParameters(Teuchos::RCP<Teuchos::ParameterList> parameterList);

    bool RandomTopologyFlag_ = false;
    bool WarmStartingFlag_ = false;
  };
}  // namespace MIRCO



// pass into MIRCO::Evaluate() just this object instead of all the params again

#endif  // SRC_INPUTPARAMETERS_H_
