#include "mirco_exportvisualization.h"

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include <cstring>  // memset
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MIRCO
{

  void ExportVisualization(const std::string& path, float gridSize, const ViewVectorInt_d activeSet,
      const std::vector<ViewMatrix_d>& otherFields, const std::vector<std::string>& otherFieldNames)
  {
    const int n = otherFields[0].extent(0);
    const int n2 = n * n;

    vtkNew<vtkImageData> img;
    img->SetDimensions(n, n, 1);
    img->SetSpacing(gridSize, gridSize, 1.0);
    img->SetOrigin(0.0, 0.0, 0.0);

    // active set
    {
      ViewVectorInt_h activeSet_h =
          Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), activeSet);

      vtkNew<vtkUnsignedCharArray> vtkArr;
      vtkArr->SetName("Active Set");
      vtkArr->SetNumberOfComponents(1);
      vtkArr->SetNumberOfTuples(n2);

      // Fast zero-fill
      auto* vtkArrArr = static_cast<unsigned char*>(vtkArr->WriteVoidPointer(0, n2));
      std::memset(vtkArrArr, 0, static_cast<size_t>(n2));

      // Unordered scatter write (duplicates don’t matter)
      for (int indA = 0; indA < activeSet_h.extent(0); ++indA)
      {
        const int a = activeSet_h(indA);
        vtkArr->SetValue(a, 1);
      }

      img->GetPointData()->AddArray(vtkArr);

      // Optional: make the first field the "active" scalar for coloring/warp defaults
      img->GetPointData()->SetActiveScalars("Active Set");
    }

    // others
    for (int f = 0; f < otherFields.size(); ++f)
    {
      ViewMatrix_h field_h =
          Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), otherFields[f]);

      vtkNew<vtkFloatArray> vtkArr;
      vtkArr->SetName(otherFieldNames[f].c_str());
      vtkArr->SetNumberOfComponents(1);
      vtkArr->SetNumberOfTuples(n2);

      // Unordered scatter write (duplicates don’t matter)
      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < n; ++j)
        {
          vtkArr->SetValue(i + n * j, field_h(i, j));
        }
      }

      img->GetPointData()->AddArray(vtkArr);
    }



    vtkNew<vtkXMLImageDataWriter> w;
    w->SetFileName((path + ".vti").c_str());
    w->SetInputData(img);

    // Binary appended data + zlib compression
    vtkNew<vtkZLibDataCompressor> z;
    w->SetCompressor(z);
    w->SetDataModeToAppended();  // binary appended payload
    w->EncodeAppendedDataOff();  // raw appended binary (not base64)

    if (!w->Write())
    {
      throw std::runtime_error("vtkXMLImageDataWriter::Write() failed");
    }
  }


}  // namespace MIRCO
