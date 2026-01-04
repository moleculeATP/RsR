#include "Recon.h"
#include <boost/filesystem.hpp>

#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

polyscope::PointCloud* ps_cloud_for_polyscope = nullptr;
polyscope::CurveNetwork* init_graph_for_polyscope = nullptr;
polyscope::CurveNetwork* edges_passed = nullptr;
std::vector<int> test_results;
std::vector<int> order_inserted;

bool SING = true;

/**
 * @brief Read the input config file and initialize
 *
 * @param config_path: Path to the config file
 * @return None
 */
void Reconstructor::read_config(fs::path config_path) {
    fs::ifstream config(config_path);
    string line;
    if (!fs::exists(config_path)) {
        throw std::runtime_error("Error opening config file!");
    }
    std::vector<string> tokens;
    while (getline(config, line)) {
        std::size_t pos = 0, found;
        while (true) {
            found = line.find(' ', pos);
            if (found != std::string::npos) {
                std::string token = line.substr(pos, found - pos);
                tokens.push_back(token);
                pos = found + 1;
            }
            else {
                found = line.find('\n', pos);
                std::string token = line.substr(pos, found - pos);
                tokens.push_back(token);
                break;
            }
        }
        string instruct = tokens[0];
        tokens[1].erase(tokens[1].find_last_not_of("\n") + 1);
        if (instruct == "root") {
            root_path = fs::path(tokens[1]);
        }
        if (instruct == "model_name") {
            if (tokens[1] == "all") {
                model_name = tokens[1];
                model_path = root_path;
            }
            else {
                model_name = tokens[1];
                model_path = root_path / tokens[1];
            }
        }
        if(instruct == "ground_truth_name") {
            fs::path gt_name = tokens[1];
            ground_truth_path = root_path / gt_name;

        }

        if (instruct == "out_root") {
            out_root_path = fs::path(tokens[1]);
        }
        if (instruct == "out_name") {
            out_name = tokens[1];
        }
        if (instruct == "isEuclidean") {
            isEuclidean = (tokens[1] == "true");
        }
        if (instruct == "isGTNormal") {
            isGTNormal = (tokens[1] == "true");
        }
        if (instruct == "isDebug") {
            isDebug = (tokens[1] == "true");
        }
        if (instruct == "isNoiseExp") {
            isNoiseExperiment = (tokens[1] == "true");
        }
        if (instruct == "k") {
            k = (std::stoi(tokens[1]));
        }
        if (instruct == "genus") {
            exp_genus = (std::stoi(tokens[1]));
        }
        if (instruct == "r") {
            r = (std::stof(tokens[1]));
        }
        if (instruct == "theta") {
            theta = (std::stof(tokens[1]));
        }
        if (instruct == "n") {
            n = (std::stoi(tokens[1]));
        }
        tokens.clear();
    }
}

/**
 * @brief Perform noise experiments by adding noise to vertex positions and normal
 *
 * @param None
 * @return None
 */
void Reconstructor::noise_experiment() {
    fs::path file_path = model_path / model_name;
    isGTNormal = false;
    std::string path_string = file_path.string();
    int start = path_string.find_last_of('\\');
    int str_len = path_string.find_last_of('.') - start;
    std::string this_model_name = path_string.substr(start + 1, str_len - 1);
    out_name = this_model_name;
    fs::path recon_root = out_root_path / "Recon" / noise_type;
    if (!fs::exists(recon_root))
        fs::create_directories(recon_root);
    /*fs::path gt_root = out_root_path / "GT" / noise_type;
    if (!fs::exists(gt_root))
        fs::create_directories(gt_root);*/

    // Load Original GT model
    io_system IO;
    {
        vector<Point> vertices;
        vector<Vector> normals;
        vector<Face> GT_faces;
        IO.read_obj(file_path, vertices, normals, GT_faces);
    }

    // Without noise
    isEuclidean = false;

    for (auto sigma : sigmas) {
        for (auto amplitude : amplitudes) {
            isFaceLoop = true;
            std::cout << "Experimenting sigma " + std::to_string(sigma) + " amplitude " + std::to_string(amplitude) + "..." << std::endl;
            std::string this_out_name = this_model_name + std::to_string(sigma) + "_" + std::to_string(amplitude);
            out_root_path = recon_root;
            out_name = this_out_name;
            reconstruct_single(noise_type, sigma, amplitude);
        }
    }
}

/**
 * @brief Reconstruct a single file
 *
 * @param noise_type: type of noise added for the noise experiments
 * @param sigma: the standard deviation of added noise
 * @param amplitude: the amplitude of added noise
 * 
 * @return None
 */
void Reconstructor::reconstruct_single(std::string noise_type, float sigma, float amplitude) {

    recon_timer.start("Initialization");
    // Init
    io_system IO;
    vector<Point> GT_vertices;
    vector<Face> GT_faces;
    std::vector<Point> in_vertices;
    std::vector<Vector> in_normals;

    // Read input
    {
        recon_timer.start("Import obj");
        fs::path file_path = root_path / model_name;
        
        std::string file_end = file_path.extension().string();
        if (file_end == ".obj") {
            IO.read_obj(file_path, in_vertices, in_normals, GT_faces);
            GT_vertices = in_vertices;
        }
        else if (file_end == ".ply"){
            IO.read_pc_ply(file_path, in_vertices, in_normals, GT_vertices);
        }
        else {
            std::cout << "Error file type" << std::endl;
            return;
        }

        if (k >= in_vertices.size())
            k = in_vertices.size() - 1;

        recon_timer.end("Import obj");
    }

    // Estimate normals & orientation & weighted smoothing
    recon_timer.start("Estimate normals");
    vector<Point> in_smoothed_v;
    {
        std::vector<int> indices(in_vertices.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(in_vertices.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(in_vertices.end(), indices.end())));
        Distance tr_dist;
        float diagonal_length;

        if (isGTNormal && in_normals.size() == 0) {
            std::cout << "No normal can be used!" << std::endl;
            isGTNormal = false;
        }

        std::vector<int> zero_normal_id;
        estimate_normal(in_vertices, kdTree, tr_dist, isGTNormal, in_normals,
            zero_normal_id, diagonal_length);

        // Fix zero normal
        if (zero_normal_id.size() > 0) {
            replace_zero_normal(zero_normal_id, kdTree, tr_dist, in_vertices, in_normals);
        }

        // Add position noise
        {
            if (amplitude != 0.) {
                float avg_edge_length = 0.0;
                for (auto face : GT_faces) {
                    for (int i = 0; i < 3; i++) {
                        Point vertex1 = in_vertices.at(i);
                        Point vertex2 = in_vertices.at((i + 1) % 3);
                        float length = norm(vertex1 - vertex2);
                        avg_edge_length += length;
                    }
                }
                avg_edge_length /= GT_faces.size() * 3;
                add_noise(noise_type, in_vertices, in_normals, sigma, amplitude, avg_edge_length);
            }
        }

        if (true) {
            std::cout << "Start first round smoothing ..." << std::endl;
            if (!isEuclidean)
                weighted_smooth(in_vertices, in_smoothed_v, in_normals, kdTree, tr_dist, diagonal_length);
            else
                in_smoothed_v = in_vertices;

            Tree kdTree1(boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.begin(), indices.begin())),
                boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.end(), indices.end())));
            Distance tr_dist1;

            estimate_normal(in_smoothed_v, kdTree1, tr_dist1, isGTNormal, in_normals,
                zero_normal_id, diagonal_length);
            if (isDebug)
                IO.export_obj(in_smoothed_v, in_normals, out_root_path / (out_name + "_smoothed_v.obj"));

            // Another round of smoothing
            if (true) {
                if (!isEuclidean) {
                    std::cout << "Start second round smoothing ..." << std::endl;
                    std::vector<Point> temp(in_smoothed_v.begin(), in_smoothed_v.end());
                    in_smoothed_v.clear();
                    weighted_smooth(temp, in_smoothed_v, in_normals, kdTree1, tr_dist1, diagonal_length);

                    Tree kdTree2(boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.begin(), indices.begin())),
                        boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.end(), indices.end())));
                    Distance tr_dist2;

                    estimate_normal(in_smoothed_v, kdTree2, tr_dist2, isGTNormal, in_normals,
                        zero_normal_id, diagonal_length);
                    if (isDebug)
                        IO.export_obj(in_smoothed_v, in_normals, out_root_path / (out_name + "_smoothed_v2.obj"));
                }
            }
        }
        else {
            in_smoothed_v = in_vertices;
        }
        
        // Add normal noise
        if(false){
            float angle = 25. / 180. * CGAL_PI;
            add_normal_noise(angle, in_normals);
        }        


    }
    recon_timer.end("Estimate normals");
    recon_timer.end("Initialization");

    recon_timer.start("algorithm");
    std::cout << "find components" << std::endl;
    // Find components
    std::vector<std::vector<Point>> component_vertices;
    std::vector<std::vector<Point>> component_smoothed_v;
    std::vector<vector<Vector>> component_normals;
    {
        // Build kdTree - CGAL
        vector<int> indices(in_smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.end(), indices.end())));
        Distance tr_dist;

        // Correct orientation

        // TODO: Seems not considerring the number of connected components when correct orientation!!!

        if (!isGTNormal) {
            s_Graph g_angle(in_vertices.size());
            s_weightMap weightmap_a = boost::get(boost::edge_weight, g_angle);
            // Init angle based graph
            for (int i = 0; i < in_vertices.size(); i++) {
                Point vertex = in_vertices[i];
                Vector this_normal = in_normals[i];

                std::vector<int> neighbors;
                std::vector<float> dists;
                kNN_search(i, vertex, kdTree, tr_dist, k, neighbors, dists);
                for (int j = 0; j < neighbors.size(); j++) {
                    if (boost::edge(i, neighbors[j], g_angle).second)
                        continue;
                    Vector neighbor_normal = in_normals[neighbors[j]];
                    boost::graph_traits< s_Graph >::edge_descriptor e;
                    bool inserted;
                    boost::tie(e, inserted) = boost::add_edge(i, neighbors[j], g_angle);
                    float angle_weight = cal_angle_based_weight(this_normal, neighbor_normal);
                    if (i == neighbors[j]) {
                        std::cout << "error" << std::endl;
                    }
                    if (angle_weight < 0)
                        std::cout << "error" << std::endl;
                    weightmap_a[e] = angle_weight;
                }
            }
            std::vector<int> component_id(in_vertices.size());
            int num = boost::connected_components(g_angle, &component_id[0]);
            //std::cout << num << std::endl;

            boost::property_map< s_Graph, boost::vertex_distance_t >::type distance_a
                = boost::get(boost::vertex_distance, g_angle);
            boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap_a
                = boost::get(boost::vertex_index, g_angle);
            std::vector<boost::graph_traits< s_Graph >::vertex_descriptor>
                p_angle(boost::num_vertices(g_angle));
            boost::prim_minimum_spanning_tree
            (g_angle, *boost::vertices(g_angle).first, &p_angle[0],
                distance_a, weightmap_a, indexmap_a,
                boost::default_dijkstra_visitor());
            g_angle.clear();

            g_angle = s_Graph(in_vertices.size());
            for (int i = 0; i < p_angle.size(); i++) {
                if (p_angle[i] == i) {
                    continue;
                }
                else {
                    boost::add_edge(i, p_angle[i], g_angle);
                }
            }
            correct_normal_orientation(g_angle, in_normals);
        }

        find_components(in_vertices, component_vertices, in_smoothed_v,
            component_smoothed_v, in_normals, component_normals, kdTree,
            tr_dist, k, isEuclidean, theta, r);

        in_vertices.clear();
        in_normals.clear();
        in_smoothed_v.clear();


        
    }

    for (int component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " + std::to_string(component_id) + " ..." << std::endl;

        isFaceLoop = true;
        std::vector<Face> faces;
        std::vector<Point> vertices = component_vertices[component_id];
        std::vector<Vector> normals = component_normals[component_id];
        std::vector<Point> smoothed_v = component_smoothed_v[component_id];


        
        // ====== POLYSCOPE VISUALIZATION ======
        {
        std::vector<std::array<double, 3>> ps_normals;
        ps_normals.reserve(normals.size());
        for (const auto& n : normals) {
            ps_normals.push_back({ n.x(), n.y(), n.z() });
        }

        std::vector<std::array<double, 3>> ps_vertices;
        ps_vertices.reserve(smoothed_v.size());
        for (const auto& p : smoothed_v) {
            ps_vertices.push_back({ p.x(), p.y(), p.z() });
        }

        ps_cloud_for_polyscope = polyscope::registerPointCloud("input point cloud", ps_vertices);
        ps_cloud_for_polyscope->setPointRadius(0.002, true);
        auto* q_normals = ps_cloud_for_polyscope->addVectorQuantity("normals", ps_normals, polyscope::VectorType::STANDARD);

        q_normals->setVectorLengthScale(0.004, true);
        q_normals->setVectorRadius(0.0005, true);
        }
        
        // ======================================





        std::string out_component_name = out_name +
            "_component_" + std::to_string(component_id);

        vector<int> indices(smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(smoothed_v.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(smoothed_v.end(), indices.end())));
        Distance tr_dist;

        recon_timer.start("Build MST");

        std::cout << "Init mst" << std::endl;

        // Initial Structure
        m_Graph mst;
        mst.graph = Graph(vertices.size());
        std::vector<m_Edge> full_edges;
        std::vector<m_Edge_length> edge_length;
        std::vector<float> connection_max_length(vertices.size(), 0.);
        std::vector<float> pre_max_length(vertices.size(), 0.);
        mst.isEuclidean = isEuclidean;
        mst.exp_genus = exp_genus;
        float max_euclidian_distance = 0.;
        {

        // ANISOTROPE SING GRAPH PART
        s_Graph g_sing;
        s_weightMap weightmap_sing = boost::get(boost::edge_weight, g_sing);
        // s_weightMap euclidian_weightmap_sing = boost::get(boost::edge_weight, g_sing);
        m_Graph sing_mst;
        float max_euclidian_distance_sing;
        std::vector<float> connection_max_length_sing(vertices.size(), 0.); // useless but testing
        std::vector<float> pre_max_length_sing(vertices.size(), 0.);

        init_sing_graph(smoothed_v, normals,
            g_sing, weightmap_sing, connection_max_length_sing, 
            25, .75, 0, exp_genus, pre_max_length_sing, max_euclidian_distance_sing, true);
        // 25 for 2000 regular
        std::cout << "Max euclidian distance in SING graph: " << max_euclidian_distance_sing << std::endl;

        // ====== POLYSCOPE VISUALIZATION : INIT SING GRAPH ======
        {
        std::vector<std::array<double,3>> ps_vertices;
        ps_vertices.reserve(smoothed_v.size());
        std::cout << smoothed_v.size() << std::endl;
        for (const auto& p : smoothed_v) {
            ps_vertices.push_back({p.x(), p.y(), p.z()});
        }

        // Convert FULL_EDGES -> polyscope edges
        std::vector<std::array<size_t,2>> ps_edges;
        ps_edges.reserve(full_edges.size());

        for (const auto& e : g_sing.m_edges) {
            ps_edges.push_back({
                static_cast<size_t>(e.m_source),
                static_cast<size_t>(e.m_target)
            });
        }

        // Register curve network
        auto* sing_graph_ =
            polyscope::registerCurveNetwork("Sing graph", ps_vertices, ps_edges);

        // Appearance
        sing_graph_->setRadius(0.0002, true);
        sing_graph_->setColor({0., 0., 0.});
        sing_graph_->setEnabled(false);
        }
            
        // ====== END POLYSCOPE VISUALIZATION ======

        // build a mst from sing edges with weightmap_sing distances
        {
            std::vector<boost::graph_traits< s_Graph >::vertex_descriptor> p_sing(boost::num_vertices(g_sing));

            boost::property_map< s_Graph, boost::vertex_distance_t >::type distance_sing
                = boost::get(boost::vertex_distance, g_sing);
            boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap_sing
                = boost::get(boost::vertex_index, g_sing);

            boost::prim_minimum_spanning_tree
            (g_sing, *boost::vertices(g_sing).first, &p_sing[0],
                distance_sing, weightmap_sing, indexmap_sing,
                boost::default_dijkstra_visitor());

            sing_mst.graph = Graph(vertices.size());
            sing_mst.isEuclidean = isEuclidean;
            sing_mst.exp_genus = exp_genus;
            build_mst(sing_mst, p_sing, isEuclidean, smoothed_v, normals);

        // ====== POLYSCOPE VISUALIZATION : sing_MST
            {
            // Convert points -> polyscope vertices
            std::vector<std::array<double,3>> ps_vertices;
            ps_vertices.reserve(smoothed_v.size());
            for (const auto& p : smoothed_v) {
                ps_vertices.push_back({p.x(), p.y(), p.z()});
            }

            // Convert FULL_EDGES -> polyscope edges
            std::vector<std::array<size_t,2>> ps_edges;
            ps_edges.reserve(boost::num_edges(sing_mst.graph));

            boost::graph_traits<Graph>::edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = boost::edges(sing_mst.graph); ei != ei_end; ++ei) {
                auto e = *ei;
                size_t u = boost::source(e, sing_mst.graph);
                size_t v = boost::target(e, sing_mst.graph);
                ps_edges.push_back({u, v});
            }

            // Register curve network
            auto* sing_mst_ =
                polyscope::registerCurveNetwork("Sing MST", ps_vertices, ps_edges);

            // Appearance
            sing_mst_->setRadius(0.001, true);
            sing_mst_->setColor({1., 0., 0.});
            sing_mst_->setEnabled(false);
            }
            // ======================================
            }

            // REGULAR SING GRAPH PART
            s_Graph g_sing_2;
            s_weightMap weightmap_sing_2 = boost::get(boost::edge_weight, g_sing);
            // s_weightMap euclidian_weightmap_sing_2 = boost::get(boost::edge_weight, g_sing_2);
            m_Graph sing_mst_2;
            float max_euclidian_distance_sing_2;
            std::vector<float> connection_max_length_sing_2(vertices.size(), 0.); // useless but testing
            std::vector<float> pre_max_length_sing_2(vertices.size(), 0.);
            
            init_sing_graph(smoothed_v, normals,
            g_sing_2, weightmap_sing_2, connection_max_length_sing_2,
            2, 1., 0, exp_genus, pre_max_length_sing_2, max_euclidian_distance_sing_2 , false);

            // ====== POLYSCOPE VISUALIZATION : INIT SING GRAPH ======
            {
            std::vector<std::array<double,3>> ps_vertices;
            ps_vertices.reserve(smoothed_v.size());
            std::cout << smoothed_v.size() << std::endl;
            for (const auto& p : smoothed_v) {
                ps_vertices.push_back({p.x(), p.y(), p.z()});
            }

            // Convert FULL_EDGES -> polyscope edges
            std::vector<std::array<size_t,2>> ps_edges;
            ps_edges.reserve(full_edges.size());

            for (const auto& e : g_sing_2.m_edges) {
                ps_edges.push_back({
                    static_cast<size_t>(e.m_source),
                    static_cast<size_t>(e.m_target)
                });
            }

            // Register curve network
            auto* sing_mst_2 =
                polyscope::registerCurveNetwork("Sing graph 2", ps_vertices, ps_edges);

            // Appearance
            sing_mst_2->setRadius(0.0002, true);
            sing_mst_2->setColor({0., 0., 0.});
            sing_mst_2->setEnabled(false);
            }
                
            // ====== END POLYSCOPE VISUALIZATION ======

            // build a mst from sing edges with weightmap_sing distances
            {
                std::vector<boost::graph_traits< s_Graph >::vertex_descriptor> p_sing(boost::num_vertices(g_sing_2));

                boost::property_map< s_Graph, boost::vertex_distance_t >::type distance_sing
                    = boost::get(boost::vertex_distance, g_sing_2);
                boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap_sing
                    = boost::get(boost::vertex_index, g_sing_2);

                boost::prim_minimum_spanning_tree
                (g_sing_2, *boost::vertices(g_sing_2).first, &p_sing[0],
                    distance_sing, weightmap_sing_2, indexmap_sing,
                    boost::default_dijkstra_visitor());

                sing_mst_2.graph = Graph(vertices.size());
                sing_mst_2.isEuclidean = isEuclidean;
                sing_mst_2.exp_genus = exp_genus;
                build_mst(sing_mst_2, p_sing, isEuclidean, smoothed_v, normals);

            // ====== POLYSCOPE VISUALIZATION : sing_MST
                {
                // Convert points -> polyscope vertices
                std::vector<std::array<double,3>> ps_vertices;
                ps_vertices.reserve(smoothed_v.size());
                for (const auto& p : smoothed_v) {
                    ps_vertices.push_back({p.x(), p.y(), p.z()});
                }

                // Convert FULL_EDGES -> polyscope edges
                std::vector<std::array<size_t,2>> ps_edges;
                ps_edges.reserve(boost::num_edges(sing_mst_2.graph));

                boost::graph_traits<Graph>::edge_iterator ei, ei_end;
                for (boost::tie(ei, ei_end) = boost::edges(sing_mst_2.graph); ei != ei_end; ++ei) {
                    auto e = *ei;
                    size_t u = boost::source(e, sing_mst_2.graph);
                    size_t v = boost::target(e, sing_mst_2.graph);
                    ps_edges.push_back({u, v});
                }

                // Register curve network
                auto* mst_2_net =
                    polyscope::registerCurveNetwork("Sing MST 2", ps_vertices, ps_edges);

                // Appearance
                mst_2_net->setRadius(0.001, true);
                mst_2_net->setColor({0., 1., 0.});
                mst_2_net->setEnabled(false);
                }
                // ======================================
            }

            // REGULAR GRAPH PART
            s_Graph g;
            s_weightMap weightmap = boost::get(boost::edge_weight, g);
            init_graph(smoothed_v, smoothed_v, normals,
                kdTree, tr_dist, k,
                g, weightmap, isEuclidean, connection_max_length,
                exp_genus, pre_max_length, theta, max_euclidian_distance);

            // ====== POLYSCOPE VISUALIZATION : INIT GRAPH ======
            {
            std::vector<std::array<double,3>> ps_vertices;
            std::vector<std::array<size_t,2>> ps_edges;
            ps_vertices.reserve(smoothed_v.size());
            std::cout << smoothed_v.size() << std::endl;
            for (const auto& p : smoothed_v) {
                ps_vertices.push_back({p.x(), p.y(), p.z()});
            }
            
            for (const auto& e : g.m_edges) {
                ps_edges.push_back({
                    static_cast<size_t>(e.m_source),
                    static_cast<size_t>(e.m_target)
                });
            }

            // Register curve network
            auto* init_graph_ =
                polyscope::registerCurveNetwork("Init graph base", ps_vertices, ps_edges);
            // Appearance
            init_graph_->setRadius(0.0002, true);
            init_graph_->setColor({0., 0., 0.});
            init_graph_->setEnabled(false);
            }
            // ======================================



            // init_sing_graph(smoothed_v, normals,
            //     g, weightmap, connection_max_length,
            //     100, 1.5, 0, exp_genus, pre_max_length);

             if(isDebug)
                IO.export_graph(g, out_root_path / (out_component_name + "_init_graph.obj"), smoothed_v, weightmap);


            // Generate MST
            std::vector<boost::graph_traits< s_Graph >::vertex_descriptor> p(boost::num_vertices(g));

            boost::property_map< s_Graph, boost::vertex_distance_t >::type distance
                = boost::get(boost::vertex_distance, g);
            boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap
                = boost::get(boost::vertex_index, g);

            boost::prim_minimum_spanning_tree
            (g, *boost::vertices(g).first, &p[0],
                distance, weightmap, indexmap,
                boost::default_dijkstra_visitor());

            build_mst(mst, p, isEuclidean, smoothed_v, normals);

        // ====== POLYSCOPE VISUALIZATION : MST ======

            std::vector<std::array<double,3>> mst_vertices;
            mst_vertices.reserve(smoothed_v.size());
            for (const auto& p : smoothed_v) {
                mst_vertices.push_back({p.x(), p.y(), p.z()});
            }

            std::vector<std::array<size_t,2>> mst_edges;
            mst_edges.reserve(boost::num_edges(mst.graph));

            boost::graph_traits<Graph>::edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = boost::edges(mst.graph); ei != ei_end; ++ei) {
                auto e = *ei;
                size_t u = boost::source(e, mst.graph);
                size_t v = boost::target(e, mst.graph);
                mst_edges.push_back({u, v});
            }

            auto* mst_net = polyscope::registerCurveNetwork("MST", mst_vertices, mst_edges);

            // Appearance settings
            mst_net->setRadius(0.001, true);
            mst_net->setColor({0., 0., 1.});
            mst_net->setEnabled(false); 

        // Store original euclidean parameters for geometry validation
        s_weightMap euclidian_dist_weightmap = weightmap;
        float euclidean_max_distance = max_euclidian_distance;
        std::vector<float> euclidean_pre_max_length = pre_max_length;
        
        if (SING){
            g = g_sing;
            mst = sing_mst;
            weightmap = weightmap_sing;  // Use SING distances for edge priority
            pre_max_length = pre_max_length_sing;  // Use SING pre_max_length for SING edge filtering
            // Use SING's actual maximum Euclidean distance for geometry validation
            euclidean_max_distance = max_euclidian_distance_sing;
            std::cout << "SING mode: using max_euclidean=" << euclidean_max_distance 
                      << " (from SING graph)" << std::endl;
        } else {
            std::cout << "Regular mode: using max_euclidean=" << euclidean_max_distance 
                      << " (from k-NN graph)" << std::endl;
        }

        // Edge arrays and sort 
            if (true) {
                int idx = 0;
                s_Graph::edge_iterator ei, ei_end;
                for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ei++) {
                    // Apply appropriate distance filtering based on mode
                    bool should_include = true;
                    if (SING) {
                        // For SING edges: filter by actual Euclidean distances to prevent spatially invalid edges
                        Point p1 = smoothed_v[(*ei).m_source];
                        Point p2 = smoothed_v[(*ei).m_target];
                        float euclidean_length = std::sqrt((p1 - p2).squared_length());
                        // Use conservative Euclidean filtering based on regular graph's pre_max_length
                        if (euclidean_length > euclidean_pre_max_length[(*ei).m_source] ||
                            euclidean_length > euclidean_pre_max_length[(*ei).m_target]) {
                            should_include = false;
                        }
                    } else {
                        // For regular edges: use original filtering
                        if (weightmap[*ei] > pre_max_length[(*ei).m_source] ||
                            weightmap[*ei] > pre_max_length[(*ei).m_target]) {
                            should_include = false;
                        }
                    }
                    
                    if (should_include) {
                        edge_length.push_back(m_Edge_length(weightmap[*ei], full_edges.size()));
                        full_edges.push_back(m_Edge((*ei).m_source, (*ei).m_target));
                        idx++;
                    }
                }
                std::sort(edge_length.begin(), edge_length.end(), edge_comparator);
                
                // Use the appropriate max distance for geometry test
                max_euclidian_distance = euclidean_max_distance;
            }
        }

        
        recon_timer.end("Build MST");

        // Export MST
        if (isDebug) {
            fs::path out_path = out_root_path / ("MST_" + out_component_name + "_C.obj");
            IO.export_graph(mst, out_path);
        }

        // ====== POLYSCOPE VISUALIZATION : INIT GRAPH ======
            {
            // Convert points -> polyscope vertices
            std::vector<std::array<double,3>> ps_vertices;
            ps_vertices.reserve(smoothed_v.size());
            for (const auto& p : smoothed_v) {
                ps_vertices.push_back({p.x(), p.y(), p.z()});
            }

            // Convert FULL_EDGES -> polyscope edges
            std::vector<std::array<size_t,2>> ps_edges;
            ps_edges.reserve(full_edges.size());

            for (const auto& e : full_edges) {
                ps_edges.push_back({
                    static_cast<size_t>(e.first),
                    static_cast<size_t>(e.second)
                });
            }

            // Register curve network
            init_graph_for_polyscope =
                polyscope::registerCurveNetwork("Init graph used", ps_vertices, ps_edges);

            // Appearance
            init_graph_for_polyscope->setRadius(0.0002, true);
            init_graph_for_polyscope->setColor({0., 0., 0.});
            init_graph_for_polyscope->setEnabled(false);
            }
            // ======================================


        






        






        // Initialize face loop label
        mst.etf.reserve(6 * vertices.size() - 11);
        init_face_loop_label(mst);

        // Betti number changes
        std::vector<int> betti_1;

        // Vanilla MST imp
        bettiNum_1 = 0;
        betti_1.push_back(bettiNum_1);
        //int inserted_edge = 0;
        if (true)
        {
            std::vector<Point> geometry_test_middle_points;
            std::vector<Vector> geometry_test_normals;

            // Custom vector for rememnering which test failed
            test_results.resize(edge_length.size(), -1);
             // -1: not tested, 0: passed, 1: topology failed, 2: geometry failed


            // Edge connection
            for (int i = 0; i < edge_length.size(); i++) {
                if (i % 100000 == 0) {
                    //std::cout << "Step " << i << " / " << edge_length.size() << std::endl;
                    showProgressBar(i / float(edge_length.size()));
                }
                unsigned int edge_idx = edge_length[i].second;
                m_Edge this_edge = full_edges[edge_idx];

                if (boost::edge(this_edge.first, this_edge.second, mst.graph).second){
                    continue; // already connected in MST
                }

                test_result_polyscope = -1;
                middle_point_geometry_test = Point(0., 0., 0.);
                normal_geometry_test = Vector(0., 0., 0.);
                bool isValid = Vanilla_check(mst, this_edge, kdTree, tr_dist, max_euclidian_distance);
                test_results[edge_idx] = test_result_polyscope;

                geometry_test_middle_points.push_back(middle_point_geometry_test);
                geometry_test_normals.push_back(normal_geometry_test);

                
                if (isValid) {
                    bool isAdded = register_face(mst, this_edge.first, this_edge.second, faces, kdTree, tr_dist, edge_length[i].first);
                    if(isAdded)
                        betti_1.push_back(bettiNum_1);
                    //Script to make FF video
                    /*inserted_edge++;
                    if (inserted_edge % int(86717 / 96) == 0) {
                        fs::path out_path = out_root_path / ("Frame_" + std::to_string(int(inserted_edge/int(86717 /96))) + ".obj");
                        IO.export_obj(vertices, mst, out_path, faces);
                        out_path = out_root_path / ("Graph_Frame_" + std::to_string(int(inserted_edge / int(86717 / 96))) + ".obj");
                        IO.export_graph(mst, out_path);
                    }*/
                }
            }
            showProgressBar(1.0);
            std::cout << std::endl;
            // std::cout << inserted_edge << std::endl;

            // Output
            if (exp_genus != 0 && isDebug) {
                fs::path out_path = out_root_path / ("Beforehandle_" + out_component_name + ".obj");
                IO.export_obj(vertices, mst, out_path, faces);
                out_path = out_root_path / ("Graph_beforehandle_" + out_component_name + ".obj");
                IO.export_graph(mst, out_path);
            }

            std::cout << "Betti 1 before handles: " << bettiNum_1 << std::endl;
            // ====== POLYSCOPE VISUALIZATION : Mesh before handles ======
            {
            // Convert vertices
            std::vector<std::array<double,3>> ps_vertices;
            ps_vertices.reserve(vertices.size());
            for (const auto& p : vertices) {
                ps_vertices.push_back({p.x(), p.y(), p.z()});
            }

            // Convert faces
            std::vector<std::array<size_t,3>> ps_faces;
            ps_faces.reserve(faces.size());
            for (const auto& f : faces) {
                ps_faces.push_back({static_cast<size_t>(f.ids[0]),
                                    static_cast<size_t>(f.ids[1]),
                                    static_cast<size_t>(f.ids[2])});
            }

            // Register the mesh in Polyscope
            auto* mesh = polyscope::registerSurfaceMesh("Graph Mesh before handles", ps_vertices, ps_faces);

            // Optional appearance
            mesh->setSurfaceColor({0.8, 0.8, 0.8});
            mesh->setEdgeWidth(0.001);
            mesh->setEnabled(false);  
            }
            // ======================================

            // ===== POLYSCOPE VISUALIZATION : Edge results =====
            {
                // Convert points -> polyscope vertices
                std::vector<std::array<double, 3>> ps_vertices;
                ps_vertices.reserve(smoothed_v.size());
                for (const auto& p : smoothed_v) {
                    ps_vertices.push_back({ p.x(), p.y(), p.z() });
                }

                // Convert FULL_EDGES -> polyscope edges with different colors according to test_results
                std::vector<std::array<size_t, 2>> ps_edges_passed;
                std::vector<std::array<size_t, 2>> ps_edges_topology_failed;
                std::vector<std::array<size_t, 2>> ps_edges_geometry_failed;
                std::vector<std::array<double, 3>> ps_middle_points;
                order_inserted.clear();
                int num_passed = 0;

                for (int i = 0; i < full_edges.size(); i++) {
                    const auto& e = full_edges[i];
                    if (test_results[i] <= 0) {
                        ps_edges_passed.push_back({
                            static_cast<size_t>(e.first),
                            static_cast<size_t>(e.second)
                            });
                        if (test_results[i] == 0){
                            order_inserted.push_back(num_passed);
                            num_passed++;
                        }
                        else {
                            order_inserted.push_back(-1); // mst edge
                        }
                    }
                    else if (test_results[i] == 1) {
                        ps_edges_topology_failed.push_back({
                            static_cast<size_t>(e.first),
                            static_cast<size_t>(e.second)
                            });
                    }
                    else if (test_results[i] == 2) {
                        ps_edges_geometry_failed.push_back({
                            static_cast<size_t>(e.first),
                            static_cast<size_t>(e.second)
                            });
                    }
                }

                // Register curve networks
                edges_passed =
                    polyscope::registerCurveNetwork("Edges Passed Before Handles", ps_vertices, ps_edges_passed);
                auto* edges_topology_failed =
                    polyscope::registerCurveNetwork("Edges Topology Failed Before Handles", ps_vertices, ps_edges_topology_failed);
                auto* edges_geometry_failed =
                    polyscope::registerCurveNetwork("Edges Geometry Failed Before Handles", ps_vertices, ps_edges_geometry_failed);

                // Appearance
                edges_passed->setRadius(0.0005, true);
                edges_passed->setColor({ 0., 1., 0. });
                edges_passed->setEnabled(false);

                edges_topology_failed->setRadius(0.0005, true);
                edges_topology_failed->setColor({ 1., 0., 0. });
                edges_topology_failed->setEnabled(false);

                edges_geometry_failed->setRadius(0.0005, true);
                edges_geometry_failed->setColor({ 0., 0., 1. });
                edges_geometry_failed->setEnabled(false);

            }
            // ======================================      
            
            // ===== POLYSCOPE VISUALIZATION : Geometry test points ======
            {
                std::vector<std::array<double, 3>> ps_test_points;
                ps_test_points.reserve(geometry_test_middle_points.size());
                for (const auto& p : geometry_test_middle_points) {
                    ps_test_points.push_back({ p.x(), p.y(), p.z() });
                }

                auto* test_points_cloud = polyscope::registerPointCloud("Geometry Test Points", ps_test_points);
                test_points_cloud->setPointRadius(0.0006, true);

                std::vector<std::array<double, 3>> ps_test_normals;
                ps_test_normals.reserve(geometry_test_normals.size());
                for (const auto& n : geometry_test_normals) {
                    ps_test_normals.push_back({ n.x(), n.y(), n.z() });
                }

                auto* q_test_normals = test_points_cloud->addVectorQuantity("Geometry Test Normals", ps_test_normals, polyscope::VectorType::STANDARD);
                q_test_normals->setVectorLengthScale(0.002, true);
                q_test_normals->setVectorRadius(0.0005, true);
                test_points_cloud->setEnabled(false);
            }
            // ======================================

        }
        // Create handles & Triangulation
        if (exp_genus != 0) {
            mst.isFinalize = true;
            std::vector<Vertex> connected_handle_root;
            connect_handle(smoothed_v, kdTree, tr_dist, max_euclidian_distance, mst, connected_handle_root, betti_1, k, isEuclidean, n);
            if (isDebug) {
                fs::path out_path = out_root_path / ("handle_" + out_component_name + ".obj");
                IO.export_edges(mst, connected_handle_root, out_path);
            }
            stop_faceloop_check();
            triangulate(faces, mst, kdTree, tr_dist, isFaceLoop, isEuclidean, connection_max_length, connected_handle_root, betti_1);
        }

        betti_1.push_back(bettiNum_1);

        // Export betti_1 : Interesting experiment with betti number.
        //IO.export_betti(betti_1, out_root_path / ("betti_" + out_component_name + ".txt"));

        // Output
        fs::path out_path = out_root_path / (out_component_name + ".obj");
        IO.export_obj(vertices, mst, out_path, faces);
        if (isDebug) {
            out_path = out_root_path / ("Graph_" + out_component_name + ".obj");
            IO.export_graph(mst, out_path, vertices);
            std::cout << "Final Betti 1: " << bettiNum_1 << std::endl;
        }

        // ====== POLYSCOPE VISUALIZATION : Mesh after handles ======

            // Convert vertices
            std::vector<std::array<double,3>> ps_vertices;
            ps_vertices.reserve(vertices.size());
            for (const auto& p : vertices) {
                ps_vertices.push_back({p.x(), p.y(), p.z()});
            }

            // Convert faces
            std::vector<std::array<size_t,3>> ps_faces;
            ps_faces.reserve(faces.size());
            for (const auto& f : faces) {
                ps_faces.push_back({static_cast<size_t>(f.ids[0]),
                                    static_cast<size_t>(f.ids[1]),
                                    static_cast<size_t>(f.ids[2])});
            }

            // Register the mesh in Polyscope
            auto* mesh = polyscope::registerSurfaceMesh("Graph Mesh after handles", ps_vertices, ps_faces);

            // Optional appearance
            mesh->setSurfaceColor({0.6, 0.8, 0.8});
            mesh->setEdgeWidth(0.001);
            mesh->setEnabled(false);  

            // ======================================
        
    }

    // ===== POLYSCOPE VISUALIZATION : ground truth
    fs::path gt_path = ground_truth_path;
    if (fs::exists(gt_path)) {
        io_system IO;
        std::vector<Point> gt_vertices;
        std::vector<Vector> gt_normals;
        std::vector<Face> gt_faces;
        std::string file_end = gt_path.extension().string();
        if (file_end == ".obj") {
            IO.read_obj(gt_path, gt_vertices, gt_normals, gt_faces);
        }
        else if (file_end == ".ply") {
            IO.read_pc_ply(gt_path, gt_vertices, gt_normals, gt_vertices);
        }

        // Convert vertices
        std::vector<std::array<double, 3>> ps_vertices;
        ps_vertices.reserve(gt_vertices.size());
        for (const auto& p : gt_vertices) {
            ps_vertices.push_back({ p.x(), p.y(), p.z() });
        }

        // Convert faces
        std::vector<std::array<size_t, 3>> ps_faces;
        ps_faces.reserve(gt_faces.size());
        for (const auto& f : gt_faces) {
            ps_faces.push_back({ static_cast<size_t>(f.ids[0]),
                                    static_cast<size_t>(f.ids[1]),
                                    static_cast<size_t>(f.ids[2]) });
        }

        // Register the mesh in Polyscope
        auto* mesh = polyscope::registerSurfaceMesh("Ground Truth Mesh", ps_vertices, ps_faces);

        // Optional appearance
        mesh->setSurfaceColor({ 0.8, 0.6, 0.6 });
        mesh->setEdgeWidth(0.001);
        mesh->setEnabled(false);
    }

    // ======================================
    
    
    recon_timer.end("algorithm");
    std::string line(40, '=');
    std::cout << line << std::endl << std::endl;
    return;
}

/**
 * @brief Entrance of reconstruction algorithm, specific timers are set
 *
 * @param None
 * @return if error happens, return -1. Otherwise return 1
 */
int Reconstructor::reconstruct() {
    // Timer whole
    recon_timer.create("Whole process");
    recon_timer.create("Initialization");
    recon_timer.create("Import obj");
    recon_timer.create("Estimate normals");
    recon_timer.create("Build MST");
    recon_timer.create("Build Rotation System");
    recon_timer.create("algorithm");

    recon_timer.start("Whole process");

    if (mode == "recon_folder") {
        traverse_and_reconstruct(model_path);
    }
    else if (mode == "single_file")
        reconstruct_single();
    else if (mode == "dtu") {
        isGTNormal = false;
        for (int i = 22; i < 129; i++) {
            std::string num = std::to_string(i);
            int n = 3;
            int precision = n - std::min<int>(n, int(num.size()));
            num = std::string(precision, '0').append(num);
            model_name = "stl" + num + "_total.ply";
            fs::path this_path(root_path / model_name);
            out_name = "stl" + num;
            if (fs::exists(this_path)) {
                std::cout << "Start Processing " + model_name + ":D" << std::endl;
                reconstruct_single();
            }
        }
    }
    else if (mode == "noise")
        noise_experiment();
    else
        return -1;
    recon_timer.end("Whole process");
    recon_timer.show();
    return 1;
}

/**
 * @brief Reconstruct a whole folder (currently only detect .obj and .ply files)
 *
 * @param dirPath: path to the directory
 * 
 * @return None
 */
void Reconstructor::traverse_and_reconstruct(const fs::path& dirPath)
{
    if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
        std::cerr << "Invalid directory: " << dirPath << std::endl;
        return;
    }

    fs::path tmp_root_path = out_root_path;
    fs::path time_recorder_path = out_root_path / "time.txt";
    fs::ofstream time_file(time_recorder_path, std::ios_base::app);
    long last_time = 0;
    for (fs::directory_iterator it(dirPath); it != fs::directory_iterator(); ++it) {
        const fs::path& filePath = it->path();
        
        if (fs::is_regular_file(filePath)) {
            if (filePath.extension() == ".obj" || filePath.extension()==".ply") {
                model_name = filePath.filename().string();
                out_root_path = tmp_root_path / filePath.stem().string();
                if (!fs::exists(out_root_path))
                    fs::create_directories(out_root_path);
                out_name = filePath.stem().string();
                reconstruct_single();
                long time = recon_timer.log("algorithm");
                time_file << out_name << ": " << time-last_time << std::endl;
                last_time = time;
            }
        }
        
        else if (fs::is_directory(filePath)) {
            traverse_and_reconstruct(filePath);
        }
    }
    time_file.close();
    return;
}

/**
 * @brief Progress bar indicating the reconstruction process
 *
 * @param progress: current progress
 *
 * @return None
 */
void Reconstructor::showProgressBar(float progress) {
    int barWidth = 70;  // Width of the progress bar

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";  // \r returns to the beginning of the line
    std::cout.flush();  // Flush the output to show the progress
}

void callback() {
    static bool is_selecting = false;
    ImGui::Checkbox("Selecting Points", &is_selecting);

    if (is_selecting) {
        ImGuiIO &io = ImGui::GetIO();
        if (io.MouseClicked[0]) {
            glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
            polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);

            if (pickResult.isHit && pickResult.structure == ps_cloud_for_polyscope) {
                int idx = ps_cloud_for_polyscope->interpretPickResult(pickResult).index;
                std::cout << "Clicked point index in v_smoothed: " << idx << std::endl;
            }else if (pickResult.isHit && pickResult.structure == init_graph_for_polyscope) {
                int idx = init_graph_for_polyscope->interpretPickResult(pickResult).index;
                
                std::vector<string> test_result_strings = {
                    "Not tested",
                    "Passed",
                    "Failed: Topology",
                    "Failed: Geometry"
                };

                std::cout
                    << "Clicked edge index in init graph: "
                    << idx << std::endl
                    << "Result of test: " << test_result_strings[test_results[idx]+1] << std::endl;
            }else if (pickResult.isHit && pickResult.structure == edges_passed){
                int idx = edges_passed->interpretPickResult(pickResult).index;
                int order = order_inserted[idx];
                std::cout << "The edge has been inserted as the " << order << "-th edge in the first phase of reconstruction." << std::endl;
            }
        }
    }
}




int main(int argc, char* argv[]){

    // Check if a parameter was provided
    std::string path;
    if (argc == 2) {
        path = argv[1]; // argv[1] is the first parameter
    }
    else if (argc > 2) {
        std::cout << "Too many arguments!" << std::endl;
        return 1;
    }
    else {
        std::cout << "Please provide the path to the config file as an argument!" << std::endl;
        return 1;
    }

    // ======== Polyscope Initialization ========
    polyscope::options::autocenterStructures = false;
    polyscope::options::autoscaleStructures = false;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();

    // ======== End of Polyscope Initialization ========

    //fs::path config_path(fs::path("../../configs/tree_config.txt"));
    fs::path config_path(path);
    Reconstructor recon(config_path);
    if (!fs::exists(recon.get_model_path())) {
        std::cout << recon.get_model_path() << " does not exist :(" << std::endl;
        return 1;
    }

    int success = recon.reconstruct();
    if (success == -1) {
        std::cout << "Mode Error!" << std::endl;
        return 1;
    }

	polyscope::state::userCallback = callback;
    polyscope::show();
    return 0;
}