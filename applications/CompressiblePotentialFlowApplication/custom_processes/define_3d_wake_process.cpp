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
                                         const Vector& rWakeNormal)
    : Process(),
      mrTrailingEdgeModelPart(rTrailingEdgeModelPart),
      mrBodyModelPart(rBodyModelPart),
      mrStlWakeModelPart(rStlWakeModelPart),
      mTolerance(Tolerance),
      mWakeNormal(rWakeNormal)
{
    KRATOS_ERROR_IF(mWakeNormal.size() != 3)
        << "The mWakeNormal should be a vector with 3 components!"
        << std::endl;
}

void Define3DWakeProcess::ExecuteInitialize()
{
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    auto& r_nodes = root_model_part.Nodes();
    VariableUtils().SetNonHistoricalVariable(WING_TIP, 0, r_nodes);

    InitializeTrailingEdgeSubModelpart();

    InitializeWakeSubModelpart();

    SetWakeAndSpanDirections();

    MarkTrailingEdgeNodes();

    ComputeLowerSurfaceNormals();

    MarkWakeElements();

    MarkKuttaElements();

    AddWakeNodes();

    CountElementNumber();

    // PrintElementIdsToFile();
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

void Define3DWakeProcess::SetWakeAndSpanDirections()
{
    const auto free_stream_velocity = mrBodyModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
    KRATOS_ERROR_IF(free_stream_velocity.size() != 3)
        << "The free stream velocity should be a vector with 3 components!"
        << std::endl;

    // Computing the norm of the free_stream_velocity vector
    const double norm = std::sqrt(inner_prod(free_stream_velocity, free_stream_velocity));

    const double eps = std::numeric_limits<double>::epsilon();
    KRATOS_ERROR_IF(norm < eps)
        << "The norm of the free stream velocity should be different than 0."
        << std::endl;

    // The wake direction is the free stream direction
    mWakeDirection = free_stream_velocity / norm;
    MathUtils<double>::CrossProduct(mSpanDirection, mWakeNormal, mWakeDirection);
}

void Define3DWakeProcess::MarkTrailingEdgeNodes()
{
    KRATOS_WATCH(mrTrailingEdgeModelPart.NumberOfNodes())
    // // chord at the root
    // const double A = 0.8104915;

    // // pi
    // const double pi = 3.14159265358979;

    // // trailing edge sweep in degrees
    // const double te_sweep = 15.69175411;
    // const double te_sweep_rad = te_sweep * pi / 180.0;

    // // angle of attack in degrees
    // const double aoa = 3.06;
    // const double aoa_rad = aoa * pi / 180.0;

    // std::vector<std::size_t> te_nodes_ordered_ids;
    // for (auto& r_node : mrBodyModelPart.Nodes()) {
    //     const double new_x = ( r_node.Y() * tan(te_sweep_rad) + A ) * cos(aoa_rad);
    //     if( r_node.X() > new_x - 0.0001){
    //         r_node.SetValue(TRAILING_EDGE, true);
    //         r_node.SetValue(KUTTA, 5.0);
    //         r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS, 0.0);
    //         r_node.SetValue(TE_ELEMENT_COUNTER, 0);
    //         te_nodes_ordered_ids.push_back(r_node.Id());
    //     }
    // }

    // std::sort(te_nodes_ordered_ids.begin(),
    //           te_nodes_ordered_ids.end());
    // mrTrailingEdgeModelPart.AddNodes(te_nodes_ordered_ids);
    // KRATOS_WATCH(mrTrailingEdgeModelPart.NumberOfNodes())

    // KRATOS_ERROR_IF(mrTrailingEdgeModelPart.NumberOfNodes() == 0) << "There are no nodes in the mrTrailingEdgeModelPart!"<< std::endl;

    /////////////////////////////////////////////////////////////////

    double max_span_position = std::numeric_limits<double>::lowest();
    double min_span_position = std::numeric_limits<double>::max();

    auto p_right_wing_tip_node = &*mrTrailingEdgeModelPart.NodesBegin();
    auto p_left_wing_tip_node = &*mrTrailingEdgeModelPart.NodesBegin();

    for (auto& r_node : mrTrailingEdgeModelPart.Nodes()) {
        r_node.SetValue(TRAILING_EDGE, true);
        r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS, 0.0);
        const auto& r_coordinates = r_node.Coordinates();
        const double distance_projection = inner_prod(r_coordinates, mSpanDirection);

        if(distance_projection > max_span_position){
            p_right_wing_tip_node = &r_node;
            max_span_position = distance_projection;
        }
        if(distance_projection < min_span_position){
            p_left_wing_tip_node = &r_node;
            min_span_position = distance_projection;
        }
    }

    mpRightWingTipNode = p_right_wing_tip_node;
    mpLeftWingTipNode = p_left_wing_tip_node;
    mpRightWingTipNode->SetValue(WING_TIP, true);
    mpLeftWingTipNode->SetValue(WING_TIP, true);
}

void Define3DWakeProcess::ComputeLowerSurfaceNormals() const
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrBodyModelPart.Conditions().size()); i++) {
        ModelPart::ConditionIterator it_cond = mrBodyModelPart.ConditionsBegin() + i;

        auto r_geometry = it_cond->GetGeometry();
        auto surface_normal = r_geometry.UnitNormal(0);
        const double projection = inner_prod(surface_normal, mWakeNormal);

        if(projection > 0.0){
            for (unsigned int i = 0; i < it_cond->GetGeometry().size(); i++) {
                r_geometry[i].SetLock();
                r_geometry[i].SetValue(NORMAL, surface_normal);
                r_geometry[i].SetValue(LOWER_SURFACE, true);
                r_geometry[i].UnSetLock();
            }
        }
        else{
            for (unsigned int i = 0; i < it_cond->GetGeometry().size(); i++) {
                r_geometry[i].SetLock();
                r_geometry[i].SetValue(UPPER_SURFACE, true);
                r_geometry[i].UnSetLock();
            }
        }
    }

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

    CalculateDistanceToSkinProcess<3> distance_calculator(root_model_part, mrStlWakeModelPart);
    distance_calculator.Execute();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(root_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = root_model_part.ElementsBegin() + i;

        // Check if the element is touching the trailing edge
        CheckIfTrailingEdgeElement(*it_elem);


        if(it_elem->Is(TO_SPLIT)){
            it_elem->SetValue(WAKE, true);
            #pragma omp critical
            {
                wake_elements_ordered_ids.push_back(it_elem->Id());
            }
            auto wake_elemental_distances_tmp = it_elem->GetValue(ELEMENTAL_DISTANCES);
            auto wake_elemental_distances = it_elem->GetValue(ELEMENTAL_DISTANCES);
            for(unsigned int j = 0; j < wake_elemental_distances.size(); j++){
                wake_elemental_distances[j] = wake_elemental_distances_tmp[j];

            }
            auto r_geometry = it_elem->GetGeometry();
            for(unsigned int j = 0; j < wake_elemental_distances.size(); j++){
                if(std::abs(wake_elemental_distances[j] < mTolerance)){
                    if(wake_elemental_distances[j] < 0.0){
                        wake_elemental_distances[j] = - mTolerance;
                    }
                    else{
                        wake_elemental_distances[j] = mTolerance;
                    }
                }
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(WAKE_DISTANCE, wake_elemental_distances[j]);
                r_geometry[j].UnSetLock();
                #pragma omp critical
                {
                    wake_nodes_ordered_ids.push_back(r_geometry[j].Id());
                }
            }
            it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
    }
    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeAndWakeElements(wake_elements_ordered_ids);
    std::sort(wake_nodes_ordered_ids.begin(),
              wake_nodes_ordered_ids.end());
    root_model_part.GetSubModelPart("wake_elements_model_part").AddNodes(wake_nodes_ordered_ids);
    KRATOS_INFO("MarkWakeElements") << "...Selecting wake elements finished..." << std::endl;
    KRATOS_CATCH("");
}

// This function checks if the element is touching the trailing edge
void Define3DWakeProcess::CheckIfTrailingEdgeElement(Element& rElement)
{
    // Loop over element nodes
    for (unsigned int i = 0; i < rElement.GetGeometry().size(); i++) {
        // Elements touching the trailing edge are trailing edge elements
        const auto& r_node = rElement.GetGeometry()[i];
        if (r_node.GetValue(TRAILING_EDGE)) {
            rElement.SetValue(TRAILING_EDGE, true);
            #pragma omp critical
            {
                mTrailingEdgeElementsOrderedIds.push_back(rElement.Id());
            }
        }
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

    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(trailing_edge_sub_model_part.Elements().size()); i++) {
        ModelPart::ElementsContainerType::ptr_iterator it_p_elem = trailing_edge_sub_model_part.Elements().ptr_begin() + i;
        auto p_elem = *it_p_elem;

        auto& r_geometry = p_elem->GetGeometry();
        // Selecting a trailing edge node
        NodeType::Pointer p_trailing_edge_node;
        unsigned int number_of_te_nodes = 0;
        for (unsigned int j = 0; j < r_geometry.size(); j++){
            const auto& r_node = r_geometry[j];
            if (r_node.GetValue(TRAILING_EDGE)){
                p_trailing_edge_node = r_geometry(j);
                number_of_te_nodes += 1;
            }
        }

        KRATOS_ERROR_IF(number_of_te_nodes < 0.5) << "Number of trailing edge nodes must be larger than 0 " << p_elem->Id()
        << " number_of_te_nodes = " << number_of_te_nodes << std::endl;

        const unsigned int number_of_non_te_nodes = 4 - number_of_te_nodes;

        Vector nodal_distances_to_te = ZeroVector(number_of_non_te_nodes);
        ComputeNodalDistancesToWakeAndLowerSurface(r_geometry, p_trailing_edge_node, nodal_distances_to_te);

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

        if(p_elem->Id()==116695){
            KRATOS_WATCH(number_of_nodes_with_negative_distance)
            KRATOS_WATCH(number_of_nodes_with_positive_distance)
        }

        // Wake structure elements (cut)
        if(number_of_nodes_with_positive_distance > 0 && number_of_nodes_with_negative_distance > 0 && p_elem->GetValue(WAKE)){
            p_elem->Set(STRUCTURE);
            BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
            unsigned int counter = 0;
            for(unsigned int j = 0; j < r_geometry.size(); j++){
                const auto& r_node = r_geometry[j];
                if(r_node.GetValue(TRAILING_EDGE)){
                    wake_elemental_distances[j] = mTolerance;
                    r_geometry[j].SetLock();
                    // Trailing edge nodes are given a positive distance
                    r_geometry[j].SetValue(WAKE_DISTANCE, mTolerance);
                    r_geometry[j].UnSetLock();
                    auto& r_number_of_neighbour_elements = r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                    #pragma omp atomic
                    r_number_of_neighbour_elements += 1;
                }
                else{
                    r_geometry[j].SetLock();
                    r_geometry[j].SetValue(WAKE_DISTANCE, nodal_distances_to_te[counter]);
                    r_geometry[j].UnSetLock();
                    wake_elemental_distances[j] = nodal_distances_to_te[counter];
                    counter += 1;
                }
            }
            p_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
        // Kutta elements (below)
        else if(number_of_nodes_with_negative_distance > number_of_non_te_nodes - 1){
            p_elem->SetValue(KUTTA, true);
            p_elem->SetValue(WAKE, false);
            p_elem->Set(TO_ERASE, true);

            // p_elem->Set(STRUCTURE);
            // p_elem->SetValue(WAKE, true);
            // #pragma omp critical
            // {
            //     wake_elements_ordered_ids.push_back(p_elem->Id());
            // }
            // BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
            // unsigned int counter = 0;
            // for(unsigned int j = 0; j < r_geometry.size(); j++){
            //     const auto& r_node = r_geometry[j];
            //     if(r_node.GetValue(TRAILING_EDGE)){
            //         wake_elemental_distances[j] = mTolerance;
            //         r_geometry[j].SetLock();
            //         r_geometry[j].SetValue(WAKE_DISTANCE, mTolerance);
            //         r_geometry[j].UnSetLock();
            //         auto& r_number_of_neighbour_elements = r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
            //         #pragma omp atomic
            //         r_number_of_neighbour_elements += 1;
            //     }
            //     else{
            //         r_geometry[j].SetLock();
            //         r_geometry[j].SetValue(WAKE_DISTANCE, nodal_distances_to_te[counter]);
            //         r_geometry[j].UnSetLock();
            //         wake_elemental_distances[j] = nodal_distances_to_te[counter];
            //         counter += 1;
            //     }
            // }
            // p_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
        }
        // Normal elements (above)
        else{
            p_elem->SetValue(WAKE, false);
            p_elem->Set(TO_ERASE, true);
        }
    }

    std::sort(wake_elements_ordered_ids.begin(),
              wake_elements_ordered_ids.end());
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");
    wake_sub_model_part.AddElements(wake_elements_ordered_ids);
    wake_sub_model_part.RemoveElements(TO_ERASE);
    KRATOS_INFO("MarkKuttaElements") << "...Selecting kutta elements finished..." << std::endl;
}

void Define3DWakeProcess::ComputeNodalDistancesToWakeAndLowerSurface(const Element::GeometryType& rGeom, NodeType::Pointer pTrailingEdgeNode, Vector& rNodalDistancesToTe) const
{
    unsigned int counter = 0;
    for (unsigned int i = 0; i < rGeom.size(); i++){
        const auto& r_node = rGeom[i];
        if (!r_node.GetValue(TRAILING_EDGE)){
            // Compute the distance vector from the trailing edge to the node
            const array_1d<double,3> distance_vector = r_node.Coordinates() - pTrailingEdgeNode->Coordinates();

            // Compute the distance in the free stream direction
            const double free_stream_direction_distance = inner_prod(distance_vector, mWakeDirection);

            double distance;
            // Nodes touching the lower surface have negative distance
            if(r_node.GetValue(LOWER_SURFACE)){
                distance = - mTolerance;
            }
            // Nodes touching the upper surface have positive distance
            else if(r_node.GetValue(UPPER_SURFACE)){
                distance = mTolerance;
            }
            else if(free_stream_direction_distance < 0.0){
                distance = inner_prod(distance_vector, pTrailingEdgeNode->GetValue(NORMAL));
                // Nodes under the wing are given a negative distance
                if(std::abs(distance) < mTolerance){
                    distance = - mTolerance;
                }
            }
            else{
                distance = inner_prod(distance_vector, mWakeNormal);
                // Nodes slightly below and above the wake are given a positive distance (wake down)
                if(std::abs(distance) < mTolerance){
                    distance = mTolerance;
                }
                // // Nodes slightly below and above the wake are given a positive distance (wake up)
                // if(std::abs(distance) < mTolerance){
                //     distance = -mTolerance;
                // }
            }

            // if(std::abs(distance) < mTolerance){
            //     distance = - mTolerance;
            // }

            rNodalDistancesToTe[counter] = distance;
            counter += 1;
        }
    }
}

void Define3DWakeProcess::AddWakeNodes() const
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

void Define3DWakeProcess::CountElementNumber() const
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

void Define3DWakeProcess::PrintElementIdsToFile() const
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
