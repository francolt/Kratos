//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Marc Nunez
//

// Project includes
#include "define_3d_wake_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "processes/calculate_distance_to_skin_process.h"

namespace Kratos {

// Constructor for Define3DWakeProcess Process
Define3DWakeProcess::Define3DWakeProcess(ModelPart& rTrailingEdgeModelPart,
                                         ModelPart& rBodyModelPart,
                                         ModelPart& rStlWakeModelPart,
                                         const double Tolerance,
                                         const Vector& rWakeNormal,
                                         const Vector& rWakeDirection,
                                         const bool SwitchWakeDirection,
                                         const bool CountElementsNumber,
                                         const bool WriteElementsIdsToFile,
                                         const bool ShedWakeFromTrailingEdge)
    : Process(),
      mrTrailingEdgeModelPart(rTrailingEdgeModelPart),
      mrBodyModelPart(rBodyModelPart),
      mrStlWakeModelPart(rStlWakeModelPart),
      mTolerance(Tolerance),
      mWakeNormal(rWakeNormal),
      mWakeDirection(rWakeDirection),
      mSwitchWakeDirection(SwitchWakeDirection),
      mCountElementsNumber(CountElementsNumber),
      mWriteElementsIdsToFile(WriteElementsIdsToFile),
      mShedWakeFromTrailingEdge(ShedWakeFromTrailingEdge)
{
    KRATOS_ERROR_IF(mWakeNormal.size() != 3)
        << "The mWakeNormal should be a vector with 3 components!"
        << std::endl;
}

void Define3DWakeProcess::ExecuteInitialize()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    auto& r_nodes = root_model_part.Nodes();
    VariableUtils().SetNonHistoricalVariable(UPPER_SURFACE, false, r_nodes);
    VariableUtils().SetNonHistoricalVariable(LOWER_SURFACE, false, r_nodes);

    InitializeTrailingEdgeSubModelpart();

    InitializeWakeSubModelpart();

    // Compute span direction as the cross product: mWakeNormal x mWakeDirection
    MathUtils<double>::CrossProduct(mSpanDirection, mWakeNormal, mWakeDirection);

    MarkTrailingEdgeNodes();

    ComputeWingLowerSurfaceNormals();

    if(mShedWakeFromTrailingEdge){
        ShedWakeSurfaceFromTheTrailingEdge();
    }

    MarkWakeElements();

    MarkKuttaElements();

    AddWakeNodesToWakeModelPart();

    if(mCountElementsNumber){
        CountElementsNumber();
    }

    if(mWriteElementsIdsToFile){
        WriteElementIdsToFile();
    }
}
// This function initializes the variables and removes all of the elements of
// the trailing edge submodelpart
void Define3DWakeProcess::InitializeTrailingEdgeSubModelpart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("trailing_edge_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // trailing_edge_sub_model_part
        ModelPart& trailing_edge_sub_model_part =
            root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

        for (auto& r_element : trailing_edge_sub_model_part.Elements()){
            r_element.SetValue(TRAILING_EDGE, false);
            r_element.SetValue(KUTTA, false);
            r_element.Reset(STRUCTURE);
            r_element.Set(TO_ERASE, true);
        }
        trailing_edge_sub_model_part.RemoveElements(TO_ERASE);
    }
    else{
        // Creating the trailing_edge_sub_model_part
        root_model_part.CreateSubModelPart("trailing_edge_elements_model_part");
    }
}

// This function initializes the variables and removes all of the elements of
// the wake submodelpart
void Define3DWakeProcess::InitializeWakeSubModelpart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("wake_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // wake_sub_model_part
        ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_elements_model_part");

        for (auto& r_element : wake_sub_model_part.Elements()){
            r_element.SetValue(WAKE, false);
            r_element.SetValue(WAKE_ELEMENTAL_DISTANCES, ZeroVector(3));
            r_element.Set(TO_ERASE, true);
        }
        wake_sub_model_part.RemoveElements(TO_ERASE);
    }
    else{
        // Creating the wake_sub_model_part
        root_model_part.CreateSubModelPart("wake_elements_model_part");
    }
}

void Define3DWakeProcess::MarkTrailingEdgeNodes() const
{
    for (auto& r_node : mrTrailingEdgeModelPart.Nodes()) {
        r_node.SetValue(TRAILING_EDGE, true);
    }
}

// This function computes the wing lower surface normals and marks the upper and
// lower surfaces. The wing lower surface normals are used later in
// ComputeNodalDistancesToWakeAndLowerSurface inside the MarkKuttaElements
// function to check whether nodes are above or below the wake
// TODO: Think a better way of doing this.
void Define3DWakeProcess::ComputeWingLowerSurfaceNormals() const
{
    // TO DISCUSS: Is it worth it to run these loops in parallel?
    // So far not much different in terms of speed has been noted.
    // Also, when ran in parallel the result is random, which I don't
    // understand why. I don't see any type of race condition here.
    // Mark upper surface
    for (auto& r_cond : mrBodyModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        const auto& surface_normal = r_geometry.UnitNormal(0);
        const double projection = inner_prod(surface_normal, mWakeNormal);

        if(!(projection > 0.0)){
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                // QUESTION: Do we need to initialize the variable UPPER_SURFACE to false?
                r_geometry[j].SetValue(UPPER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }

    // Mark lower surface
    for (auto& r_cond : mrBodyModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        const auto& surface_normal = r_geometry.UnitNormal(0);
        const double projection = inner_prod(surface_normal, mWakeNormal);

        if(projection > 0.0){
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(NORMAL, surface_normal);
                r_geometry[j].SetValue(LOWER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }
}

void Define3DWakeProcess::ShedWakeSurfaceFromTheTrailingEdge() const
{
    // ModelPart& root_model_part = mrTrailingEdgeModelPart.GetRootModelPart();
    // root_model_part.CreateSubModelPart("wake_model_part");
    // ModelPart& wake_surface_model_part = root_model_part.GetSubModelPart("wake_model_part");
    Properties::Pointer pElemProp = mrStlWakeModelPart.CreateNewProperties(0);
    const double shedded_distance = 12.5;
    const double element_size = 0.2;
    const double number_of_elements = shedded_distance / element_size;
    const unsigned int number_of_elements_in_wake_direction = int(number_of_elements);
    const double z = 1e-9;
    // KRATOS_WATCH(shedded_distance)
    KRATOS_WATCH(element_size)
    KRATOS_WATCH(number_of_elements)
    KRATOS_WATCH(number_of_elements_in_wake_direction)
    KRATOS_WATCH(z)
    IndexType node_index = 0;
    IndexType element_index = 0;

    array_1d<double,3> coordinates1 = ZeroVector(3);
    array_1d<double,3> coordinates2 = ZeroVector(3);
    array_1d<double,3> coordinates3 = ZeroVector(3);
    array_1d<double,3> coordinates4 = ZeroVector(3);

    for (auto& r_cond : mrTrailingEdgeModelPart.Conditions()) {
        const auto& r_geometry = r_cond.GetGeometry();
        coordinates1 = r_geometry[0].Coordinates();
        coordinates2 = r_geometry[1].Coordinates();
        coordinates3 = coordinates1 + element_size * mWakeDirection;
        coordinates4 = coordinates2 + element_size * mWakeDirection;

        const auto& p_node1 = mrStlWakeModelPart.CreateNewNode(
            ++node_index, coordinates1[0], coordinates1[1], coordinates1[2] + z);
        const auto& p_node2 = mrStlWakeModelPart.CreateNewNode(
            ++node_index, coordinates2[0], coordinates2[1], coordinates2[2] + z);
        const auto& p_node3 = mrStlWakeModelPart.CreateNewNode(
            ++node_index, coordinates3[0], coordinates3[1], coordinates3[2] + z);
        const auto& p_node4 = mrStlWakeModelPart.CreateNewNode(
            ++node_index, coordinates4[0], coordinates4[1], coordinates4[2] + z);

        const auto& side1 = coordinates2 - coordinates1;
        const auto& side2 = coordinates3 - coordinates1;
        array_1d<double, 3> face_normal = ZeroVector(3);
        MathUtils<double>::CrossProduct(face_normal, side1, side2);
        const double normal_projection = inner_prod(face_normal, mWakeNormal);

        if(normal_projection > 0.0){
            const std::vector<ModelPart::IndexType> elemNodes1{p_node1->Id(), p_node2->Id(), p_node3->Id()};
            const std::vector<ModelPart::IndexType> elemNodes2{p_node1->Id(), p_node4->Id(), p_node3->Id()};
            mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes1, pElemProp);
            mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes2, pElemProp);
        }
        else{
            const std::vector<ModelPart::IndexType> elemNodes1{p_node1->Id(), p_node3->Id(), p_node2->Id()};
            const std::vector<ModelPart::IndexType> elemNodes2{p_node1->Id(), p_node3->Id(), p_node4->Id()};
            mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes1, pElemProp);
            mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes2, pElemProp);
        }

        for (unsigned int j = 0; j < number_of_elements_in_wake_direction; j++){
            coordinates1 = coordinates3;
            coordinates2 = coordinates4;
            coordinates3 = coordinates1 + element_size * mWakeDirection;
            coordinates4 = coordinates2 + element_size * mWakeDirection;

            const auto& p_node1 = mrStlWakeModelPart.CreateNewNode(
                ++node_index, coordinates1[0], coordinates1[1], coordinates1[2] + z);
            const auto& p_node2 = mrStlWakeModelPart.CreateNewNode(
                ++node_index, coordinates2[0], coordinates2[1], coordinates2[2] + z);
            const auto& p_node3 = mrStlWakeModelPart.CreateNewNode(
                ++node_index, coordinates3[0], coordinates3[1], coordinates3[2] + z);
            const auto& p_node4 = mrStlWakeModelPart.CreateNewNode(
                ++node_index, coordinates4[0], coordinates4[1], coordinates4[2] + z);

            const auto& side1 = coordinates2 - coordinates1;
            const auto& side2 = coordinates3 - coordinates1;
            array_1d<double, 3> face_normal = ZeroVector(3);
            MathUtils<double>::CrossProduct(face_normal, side1, side2);
            const double normal_projection = inner_prod(face_normal, mWakeNormal);

            if(normal_projection > 0.0){
                const std::vector<ModelPart::IndexType> elemNodes1{p_node1->Id(), p_node2->Id(), p_node3->Id()};
                const std::vector<ModelPart::IndexType> elemNodes2{p_node1->Id(), p_node4->Id(), p_node3->Id()};
                mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes1, pElemProp);
                mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes2, pElemProp);
            }
            else{
                const std::vector<ModelPart::IndexType> elemNodes1{p_node1->Id(), p_node3->Id(), p_node2->Id()};
                const std::vector<ModelPart::IndexType> elemNodes2{p_node1->Id(), p_node3->Id(), p_node4->Id()};
                mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes1, pElemProp);
                mrStlWakeModelPart.CreateNewElement("Element3D3N", ++element_index, elemNodes2, pElemProp);
            }
        }
    }
    KRATOS_WATCH(mrStlWakeModelPart.NumberOfNodes())
    KRATOS_WATCH(mrStlWakeModelPart.NumberOfElements())
}

// This function checks which elements are cut by the wake and marks them as
// wake elements
void Define3DWakeProcess::MarkWakeElements()
{
    KRATOS_TRY;
    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements..." << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    std::vector<std::size_t> wake_elements_ordered_ids;
    std::vector<std::size_t> wake_nodes_ordered_ids;

    // TODO: substitute with CalculateDiscontinuousDistanceToSkinProcess
    const auto start_time(std::chrono::steady_clock::now());
    CalculateDistanceToSkinProcess<3> distance_calculator(root_model_part, mrStlWakeModelPart);
    distance_calculator.Execute();
    const auto end_time = std::chrono::duration_cast<std::chrono::duration<double>>(
                              std::chrono::steady_clock::now() - start_time)
                              .count();
    KRATOS_INFO("MarkWakeElements") << " distance_calculator took " << end_time
                                    << " [sec]" << std::endl;

    // This variable allows to inverse the distances computed with the distance
    // process, it is useful if the user makes a mistake and defines the wake
    // surface with the normal pointing downwards.
    double wake_normal_switching_factor = 1.0;
    if(mSwitchWakeDirection){
        KRATOS_INFO("MarkWakeElements") << " Switching wake element distances!" << std::endl;
        wake_normal_switching_factor = -1.0;
    }
    // std::ofstream outfile_wake;
    // outfile_wake.open("wake_elements_id.txt");

    block_for_each(root_model_part.Elements(), [&](Element& rElement)
    {
        // Check if the element is touching the trailing edge
        auto& r_geometry = rElement.GetGeometry();
        CheckIfTrailingEdgeElement(rElement, r_geometry);

        // Mark wake elements, save their ids, save the elemental distances in
        // the element and in the nodes, and save wake nodes ids
        if (rElement.Is(TO_SPLIT))
        {
            // Mark wake elements
            rElement.SetValue(WAKE, true);
            // Save wake elements ids
            #pragma omp critical
            {
                // std::ofstream outfile_wake;
                // outfile_wake.open("wake_elements_id.txt", std::ios_base::app);
                // outfile_wake << rElement.Id();
                // outfile_wake << "\n";
                wake_elements_ordered_ids.push_back(rElement.Id());
            }
            // Save elemental distances in the element
            array_1d<double,4>  wake_elemental_distances = ZeroVector(4);
            wake_elemental_distances = wake_normal_switching_factor * rElement.GetValue(ELEMENTAL_DISTANCES);
            rElement.SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
            // Save elemental distances in the nodes
            for (unsigned int j = 0; j < wake_elemental_distances.size(); j++)
            {
                if (std::abs(wake_elemental_distances[j] < mTolerance))
                {
                    if (wake_elemental_distances[j] < 0.0)
                    {
                        wake_elemental_distances[j] = - mTolerance;
                    }
                    else{
                        wake_elemental_distances[j] = mTolerance;
                    }
                }
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(WAKE_DISTANCE, wake_elemental_distances[j]);
                r_geometry[j].UnSetLock();
                // Save nodes ids
                #pragma omp critical
                {
                    wake_nodes_ordered_ids.push_back(r_geometry[j].Id());
                }
            }

        }
    });

    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeAndWakeElements(wake_elements_ordered_ids);
    std::sort(wake_nodes_ordered_ids.begin(),
              wake_nodes_ordered_ids.end());
    root_model_part.GetSubModelPart("wake_elements_model_part").AddNodes(wake_nodes_ordered_ids);
    KRATOS_WATCH(root_model_part.GetSubModelPart("wake_elements_model_part").NumberOfNodes());
    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements finished..." << std::endl;
    KRATOS_CATCH("");
}

// This function checks if the element is touching the trailing edge
void Define3DWakeProcess::CheckIfTrailingEdgeElement(Element& rElement, Geometry<NodeType>& rGeometry)
{
    // Loop over element nodes
    for (unsigned int i = 0; i < rGeometry.size(); i++) {
        // Elements touching the trailing edge are trailing edge elements
        rGeometry[i].SetLock();
        if (rGeometry[i].GetValue(TRAILING_EDGE)) {
            rElement.SetValue(TRAILING_EDGE, true);
            #pragma omp critical
            {
                mTrailingEdgeElementsOrderedIds.push_back(rElement.Id());
            }
        }
        rGeometry[i].UnSetLock();
    }
}

// This function adds the trailing edge elements in the
// trailing_edge_sub_model_part
void Define3DWakeProcess::AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds)
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    std::sort(rWakeElementsOrderedIds.begin(),
              rWakeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("wake_elements_model_part").AddElements(rWakeElementsOrderedIds);

    std::sort(mTrailingEdgeElementsOrderedIds.begin(),
              mTrailingEdgeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("trailing_edge_elements_model_part").AddElements(mTrailingEdgeElementsOrderedIds);
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
void Define3DWakeProcess::MarkKuttaElements()
{
    KRATOS_INFO("MarkKuttaElements") << "...Selecting kutta elements..." << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    std::vector<std::size_t> wake_elements_ordered_ids;

    // TO DISCUSS: Is it worth it to run this loop in parallel?
    // So far not much different in terms of speed has been noted.
    block_for_each(trailing_edge_sub_model_part.Elements(), [&](Element& rElement)
    {
        auto& r_geometry = rElement.GetGeometry();
        // Selecting one trailing edge node from the element
        // (this gets the last one according to the order within r_geometry)
        NodeType::Pointer p_trailing_edge_node;
        unsigned int number_of_te_nodes = 0;
        for (unsigned int j = 0; j < r_geometry.size(); j++){
            const auto& r_node = r_geometry[j];
            if (r_node.GetValue(TRAILING_EDGE)){
                p_trailing_edge_node = r_geometry(j);
                number_of_te_nodes += 1;
            }
        }

        KRATOS_ERROR_IF(number_of_te_nodes < 1)
            << "Number of trailing edge nodes must be 1 or larger. Element Id: "
            << rElement.Id() << " number_of_te_nodes = " << number_of_te_nodes
            << std::endl;

        const unsigned int number_of_non_te_nodes = 4 - number_of_te_nodes;

        Vector nodal_distances_to_te = ZeroVector(number_of_non_te_nodes);
        ComputeNodalDistancesToWakeOrWingLowerSurface(r_geometry, p_trailing_edge_node, nodal_distances_to_te);

        unsigned int number_of_nodes_with_negative_distance = 0;
        unsigned int number_of_nodes_with_positive_distance = 0;

        for(unsigned int j = 0; j < nodal_distances_to_te.size(); j++){
            if(nodal_distances_to_te[j] < 0.0){
                number_of_nodes_with_negative_distance += 1;
            }
            else{
                number_of_nodes_with_positive_distance +=1;
            }
        }

        // Wake structure elements (cut)
        if(number_of_nodes_with_positive_distance > 0 && number_of_nodes_with_negative_distance > 0 && rElement.GetValue(WAKE)){
            rElement.Set(STRUCTURE);
            BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
            unsigned int counter = 0;
            for(unsigned int j = 0; j < r_geometry.size(); j++){
                const auto& r_node = r_geometry[j];
                if(r_node.GetValue(TRAILING_EDGE)){
                    // Trailing edge nodes are given a positive distance
                    wake_elemental_distances[j] = mTolerance;
                    r_geometry[j].SetLock();
                    r_geometry[j].SetValue(WAKE_DISTANCE, mTolerance);
                    r_geometry[j].UnSetLock();
                }
                else{
                    // Assigning computed distance
                    r_geometry[j].SetLock();
                    r_geometry[j].SetValue(WAKE_DISTANCE, nodal_distances_to_te[counter]);
                    r_geometry[j].UnSetLock();
                    wake_elemental_distances[j] = nodal_distances_to_te[counter];
                    counter += 1;
                }
            }
            rElement.SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
        // Kutta elements (below). Kutta elements have all non_te_nodes with negative distance:
        // 1 te node  -> 3 non te nodes with negative distance
        // 2 te nodes -> 2 non te nodes with negative distance
        else if(number_of_nodes_with_negative_distance > number_of_non_te_nodes - 1){
            rElement.SetValue(KUTTA, true);
            rElement.SetValue(WAKE, false);
            rElement.Set(TO_ERASE, true);
        }
        // Normal elements (above). Normal elements have all nodes with positive distance.
        else{
            rElement.SetValue(WAKE, false);
            rElement.Set(TO_ERASE, true);
        }
    });

    std::sort(wake_elements_ordered_ids.begin(),
              wake_elements_ordered_ids.end());
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");
    wake_sub_model_part.AddElements(wake_elements_ordered_ids);
    wake_sub_model_part.RemoveElements(TO_ERASE);
    KRATOS_INFO("MarkKuttaElements") << "...Selecting kutta elements finished..." << std::endl;
}

// This function computes the distances from the current element non trailing
// edge nodes to the wake or to the wing lower surface, depending whether the
// node is behind or in front of the trailing edge according to the wake
// direction (free stream velocity direction).
void Define3DWakeProcess::ComputeNodalDistancesToWakeOrWingLowerSurface(const Element::GeometryType& rGeom, NodeType::Pointer pTrailingEdgeNode, Vector& rNodalDistancesToTe) const
{
    unsigned int counter = 0;
    for (unsigned int i = 0; i < rGeom.size(); i++){
        const auto& r_node = rGeom[i];
        if (!r_node.GetValue(TRAILING_EDGE)){
            // Compute the distance vector from the selected trailing edge node to the current node
            const array_1d<double,3> distance_vector = r_node.Coordinates() - pTrailingEdgeNode->Coordinates();

            // Compute the distance in the free stream direction
            const double free_stream_direction_distance = inner_prod(distance_vector, mWakeDirection);

            double distance;
            // Nodes in the lower surface are directly assigned a negative distance
            if(r_node.GetValue(LOWER_SURFACE)){
                distance = - mTolerance;
            }
            // Nodes in the upper surface are directly assigned a positive distance
            else if(r_node.GetValue(UPPER_SURFACE)){
                distance = mTolerance;
            }
            //  For the nodes in front of the selected trailing edge node, the
            //  distance is computed according to the lower surface normal
            else if(free_stream_direction_distance < 0.0){
                distance = inner_prod(distance_vector, pTrailingEdgeNode->GetValue(NORMAL));
                // Nodes under the wing are given a negative distance
                if(std::abs(distance) < mTolerance){
                    distance = - mTolerance;
                }
            }
            // For the nodes behind of the selected trailing edge node, the
            // distance is computed according to the wake normal
            else
            {
                distance = inner_prod(distance_vector, mWakeNormal);
                // Nodes slightly below and above the wake are given a positive distance (wake down)
                if(std::abs(distance) < mTolerance){
                    distance = mTolerance;
                }
            }
            rNodalDistancesToTe[counter] = distance;
            counter += 1;
        }
    }
}

void Define3DWakeProcess::AddWakeNodesToWakeModelPart() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_elements_model_part");

    std::vector<std::size_t> wake_nodes_ordered_ids;
    for (auto& r_element : wake_sub_model_part.Elements()){
        for (unsigned int i = 0; i < r_element.GetGeometry().size(); i++){
            r_element.GetGeometry()[i].SetValue(WAKE, true);
            wake_nodes_ordered_ids.push_back(r_element.GetGeometry()[i].Id());
        }
    }

    std::sort(wake_nodes_ordered_ids.begin(),
              wake_nodes_ordered_ids.end());
    wake_sub_model_part.AddNodes(wake_nodes_ordered_ids);
}

// This function counts the number of elements of each type. Useful for
// debugging purposes.
void Define3DWakeProcess::CountElementsNumber() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    // Initialize counters
    unsigned int normal_elements_counter = 0;
    unsigned int wake_elements_counter = 0;
    unsigned int kutta_elements_counter = 0;
    unsigned int structure_elements_counter = 0;
    for (auto& r_element : trailing_edge_sub_model_part.Elements()){
        if(!r_element.GetValue(WAKE)){
            if(r_element.GetValue(KUTTA)){
                kutta_elements_counter += 1;
            }
            else{
                normal_elements_counter += 1;
            }
        }
        else{
            wake_elements_counter += 1;
            if(r_element.Is(STRUCTURE)){
                structure_elements_counter += 1;
            }
        }
    }

    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");
    const unsigned int all_wake_elements_counter = wake_sub_model_part.NumberOfElements();
    KRATOS_WATCH(normal_elements_counter)
    KRATOS_WATCH(kutta_elements_counter)
    KRATOS_WATCH(wake_elements_counter)
    KRATOS_WATCH(structure_elements_counter)
    KRATOS_WATCH(all_wake_elements_counter)
}

// This function prints the elements ids into a file. Useful for debugging
// purposes.
void Define3DWakeProcess::WriteElementIdsToFile() const
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    std::ofstream outfile;
    outfile.open("normal_elements_id.txt");
    std::ofstream outfile_wake;
    outfile_wake.open("wake_elements_id.txt");
    std::ofstream outfile_structure;
    outfile_structure.open("structure_elements_id.txt");
    std::ofstream outfile_kutta;
    outfile_kutta.open("kutta_elements_id.txt");
        for (auto& r_element : trailing_edge_sub_model_part.Elements()){
        if(!r_element.GetValue(WAKE)){
            if(r_element.GetValue(KUTTA)){
                std::ofstream outfile_kutta;
                outfile_kutta.open("kutta_elements_id.txt", std::ios_base::app);
                outfile_kutta << r_element.Id();
                outfile_kutta << "\n";
            }
            else{
                std::ofstream outfile;
                outfile.open("normal_elements_id.txt", std::ios_base::app);
                outfile << r_element.Id();
                outfile << "\n";
            }
        }
        else{
            std::ofstream outfile_wake;
            outfile_wake.open("wake_elements_id.txt", std::ios_base::app);
            outfile_wake << r_element.Id();
            outfile_wake << "\n";
            if(r_element.Is(STRUCTURE)){
                std::ofstream outfile_structure;
                outfile_structure.open("structure_elements_id.txt", std::ios_base::app);
                outfile_structure << r_element.Id();
                outfile_structure << "\n";
            }
        }
    }

    // Loop over all wake elements
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");

    std::ofstream outfile_all_wake;
    outfile_wake.open("all_wake_elements_id.txt");
    for (auto& r_element : wake_sub_model_part.Elements()){
        std::ofstream outfile_all_wake;
        outfile_all_wake.open("all_wake_elements_id.txt", std::ios_base::app);
        outfile_all_wake << r_element.Id();
        outfile_all_wake << "\n";

    }
}
} // namespace Kratos.
