"""Interactive Streamlit dashboard for BioPilot."""

import io
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import pandas as pd
import numpy as np

import streamlit as st

from biopilot.src.fetcher import DatasetFetcher
from biopilot.src.annotation import AnnotationDB
from biopilot.src.analyzer import DataAnalyzer
from biopilot.src.reproducibility import ReproducibilityLogger

st.set_page_config(
    page_title="BioPilot Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
    }
    .subheader {
        font-size: 1.5rem;
        color: #2c3e50;
    }
    .metric-card {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 15px;
        margin: 5px 0;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_resource
def get_fetcher():
    return DatasetFetcher(data_dir="biopilot/data")


@st.cache_resource
def get_db():
    return AnnotationDB(db_path="biopilot/data/metadata/annotations.db")


@st.cache_resource
def get_analyzer():
    return DataAnalyzer(results_dir="biopilot/results", logs_dir="biopilot/logs")


@st.cache_resource
def get_logger():
    return ReproducibilityLogger(logs_dir="biopilot/logs", env_dir="biopilot/env")


def render_sidebar():
    st.sidebar.markdown("# 🧬 BioPilot")
    st.sidebar.markdown("*Bioinformatics Assistant*")
    st.sidebar.markdown("---")
    
    page = st.sidebar.radio(
        "Navigation",
        ["🏠 Home", "🔍 Data Search", "📊 Sample Manager", "📈 Analysis", "📋 Pipelines", "⚙️ Settings"],
        label_visibility="collapsed",
    )
    
    st.sidebar.markdown("---")
    
    db = get_db()
    stats = db.get_statistics()
    
    st.sidebar.markdown("### Quick Stats")
    st.sidebar.metric("Total Samples", stats.get("total_samples", 0))
    st.sidebar.metric("Species", stats.get("species_count", 0))
    
    return page


def render_home():
    st.markdown('<h1 class="main-header">🧬 BioPilot Dashboard</h1>', unsafe_allow_html=True)
    st.markdown("### A lightweight, reproducible bioinformatics assistant")
    
    st.markdown("---")
    
    col1, col2, col3, col4 = st.columns(4)
    
    db = get_db()
    stats = db.get_statistics()
    
    with col1:
        st.metric("📊 Total Samples", stats.get("total_samples", 0))
    with col2:
        st.metric("🦠 Species", stats.get("species_count", 0))
    with col3:
        st.metric("🧪 Tissue Types", len(stats.get("by_tissue", {})))
    with col4:
        st.metric("🔬 Seq Types", len(stats.get("by_sequencing_type", {})))
    
    st.markdown("---")
    
    col_left, col_right = st.columns(2)
    
    with col_left:
        st.markdown("### Samples by Species")
        by_species = stats.get("by_species", {})
        if by_species:
            df_species = pd.DataFrame(
                list(by_species.items()),
                columns=["Species", "Count"]
            ).sort_values("Count", ascending=False)
            st.bar_chart(df_species.set_index("Species"))
        else:
            st.info("No samples in database yet.")
    
    with col_right:
        st.markdown("### Samples by Sequencing Type")
        by_seq = stats.get("by_sequencing_type", {})
        if by_seq:
            df_seq = pd.DataFrame(
                list(by_seq.items()),
                columns=["Type", "Count"]
            )
            st.bar_chart(df_seq.set_index("Type"))
        else:
            st.info("No samples in database yet.")
    
    st.markdown("---")
    
    st.markdown("### Recent Activity")
    logger = get_logger()
    recent_commands = logger.commands[-5:] if logger.commands else []
    
    if recent_commands:
        for cmd in reversed(recent_commands):
            st.text(f"[{cmd.timestamp}] {cmd.command}")
    else:
        st.info("No recent activity.")


def render_search():
    st.markdown('<h1 class="subheader">🔍 Data Search</h1>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["GEO", "SRA", "ENA", "Downloaded"])
    
    fetcher = get_fetcher()
    
    with tab1:
        st.markdown("### Search GEO Database")
        
        geo_query = st.text_input("Search query", placeholder="e.g., RNA-seq human cancer")
        max_results = st.slider("Max results", 5, 50, 10)
        
        if st.button("Search GEO", type="primary"):
            if geo_query:
                with st.spinner("Searching..."):
                    results = fetcher.search_geo(geo_query, max_results)
                
                if results:
                    st.success(f"Found {len(results)} datasets")
                    
                    for ds in results:
                        with st.expander(f"**{ds.accession}**: {ds.title[:60]}..."):
                            col1, col2 = st.columns(2)
                            with col1:
                                st.markdown(f"**Organism:** {ds.organism}")
                                st.markdown(f"**Platform:** {ds.platform}")
                            with col2:
                                st.markdown(f"**Samples:** {ds.sample_count}")
                                st.markdown(f"**Type:** {ds.study_type}")
                            
                            if st.button(f"Download {ds.accession}", key=f"geo_dl_{ds.accession}"):
                                metadata_file = fetcher.save_metadata(ds)
                                st.success(f"Saved metadata to {metadata_file}")
                else:
                    st.warning("No results found.")
    
    with tab2:
        st.markdown("### Search SRA Database")
        
        sra_query = st.text_input("SRA query", placeholder="e.g., RNA-seq")
        max_sra = st.slider("Max results", 5, 50, 10, key="sra_max")
        
        if st.button("Search SRA", type="primary"):
            if sra_query:
                with st.spinner("Searching..."):
                    results = fetcher.search_sra(sra_query, max_sra)
                
                if results:
                    st.success(f"Found {len(results)} runs")
                    
                    df = pd.DataFrame([
                        {
                            "Accession": r.accession,
                            "Organism": r.organism,
                            "Platform": r.platform,
                            "Type": r.study_type,
                        }
                        for r in results
                    ])
                    st.dataframe(df, use_container_width=True)
    
    with tab3:
        st.markdown("### Search ENA Database")
        
        ena_query = st.text_input("ENA query", placeholder="e.g., Homo sapiens")
        max_ena = st.slider("Max results", 5, 50, 10, key="ena_max")
        
        if st.button("Search ENA", type="primary"):
            if ena_query:
                with st.spinner("Searching..."):
                    results = fetcher.search_ena(ena_query, max_ena)
                
                if results:
                    st.success(f"Found {len(results)} runs")
                    
                    df = pd.DataFrame([
                        {
                            "Accession": r.accession,
                            "Organism": r.organism,
                            "Platform": r.platform,
                            "Type": r.study_type,
                        }
                        for r in results
                    ])
                    st.dataframe(df, use_container_width=True)
    
    with tab4:
        st.markdown("### Downloaded Datasets")
        
        downloaded = fetcher.list_downloaded()
        
        if downloaded:
            df = pd.DataFrame(downloaded)
            st.dataframe(df, use_container_width=True)
        else:
            st.info("No downloaded datasets yet.")


def render_samples():
    st.markdown('<h1 class="subheader">📊 Sample Manager</h1>', unsafe_allow_html=True)
    
    db = get_db()
    
    tab1, tab2, tab3 = st.tabs(["View Samples", "Add Sample", "Search & Filter"])
    
    with tab1:
        samples = db.query_samples()
        
        if samples:
            df = pd.DataFrame([
                {
                    "Sample ID": s.sample_id,
                    "Accession": s.accession,
                    "Species": s.species,
                    "Tissue": s.tissue_type,
                    "Seq Type": s.sequencing_type,
                    "Condition": s.condition,
                }
                for s in samples
            ])
            
            st.dataframe(df, use_container_width=True)
        else:
            st.info("No samples in database.")
    
    with tab2:
        st.markdown("### Add New Sample")
        
        with st.form("add_sample_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                sample_id = st.text_input("Sample ID *")
                accession = st.text_input("Accession")
                species = st.text_input("Species *")
                tissue = st.text_input("Tissue Type")
            
            with col2:
                seq_type = st.selectbox(
                    "Sequencing Type",
                    ["RNA-seq", "ChIP-seq", "ATAC-seq", "WGS", "WES", "Other"]
                )
                platform = st.selectbox(
                    "Platform",
                    ["Illumina", "PacBio", "Nanopore", "Other"]
                )
                condition = st.text_input("Condition")
                replicate = st.number_input("Replicate", min_value=0, value=0)
            
            submitted = st.form_submit_button("Add Sample", type="primary")
            
            if submitted and sample_id and species:
                from biopilot.src.annotation.database import Sample
                sample = Sample(
                    sample_id=sample_id,
                    accession=accession or "",
                    species=species,
                    tissue_type=tissue or "",
                    sequencing_type=seq_type,
                    platform=platform,
                    condition=condition or "",
                    replicate=replicate if replicate > 0 else None,
                )
                db.add_sample(sample)
                st.success(f"Added sample: {sample_id}")
                st.rerun()
    
    with tab3:
        st.markdown("### Search & Filter")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            filter_species = st.selectbox(
                "Species",
                ["All"] + db.list_species()
            )
        with col2:
            filter_tissue = st.selectbox(
                "Tissue Type",
                ["All"] + db.list_tissue_types()
            )
        with col3:
            filter_seq = st.selectbox(
                "Seq Type",
                ["All"] + db.list_sequencing_types()
            )
        
        filters = {}
        if filter_species != "All":
            filters["species"] = filter_species
        if filter_tissue != "All":
            filters["tissue_type"] = filter_tissue
        if filter_seq != "All":
            filters["sequencing_type"] = filter_seq
        
        samples = db.query_samples(**filters)
        
        if samples:
            df = pd.DataFrame([
                {
                    "Sample ID": s.sample_id,
                    "Species": s.species,
                    "Tissue": s.tissue_type,
                    "Condition": s.condition,
                }
                for s in samples
            ])
            st.dataframe(df, use_container_width=True)
        else:
            st.info("No matching samples.")


def render_analysis():
    st.markdown('<h1 class="subheader">📈 Data Analysis</h1>', unsafe_allow_html=True)
    
    analyzer = get_analyzer()
    
    tab1, tab2, tab3, tab4 = st.tabs(["Load Data", "Normalization", "PCA", "Differential Expression"])
    
    with tab1:
        st.markdown("### Load Expression Data")
        
        uploaded_file = st.file_uploader(
            "Upload expression matrix (CSV/TSV)",
            type=["csv", "tsv", "txt"]
        )
        
        if uploaded_file:
            sep = "\t" if uploaded_file.name.endswith((".tsv", ".txt")) else ","
            
            try:
                df = pd.read_csv(uploaded_file, sep=sep, index_col=0)
                st.success(f"Loaded: {df.shape[0]} genes x {df.shape[1]} samples")
                
                st.session_state["expression_data"] = df
                
                st.markdown("#### Preview")
                st.dataframe(df.head(10), use_container_width=True)
                
                st.markdown("#### Summary Statistics")
                st.write(df.describe())
            except Exception as e:
                st.error(f"Error loading file: {e}")
    
    with tab2:
        st.markdown("### Normalization")
        
        if "expression_data" not in st.session_state:
            st.warning("Please load data first.")
        else:
            data = st.session_state["expression_data"]
            
            method = st.selectbox(
                "Normalization Method",
                ["cpm", "rpkm", "tpm", "zscore", "quantile", "minmax", "robust"]
            )
            log_transform = st.checkbox("Log2 transform", value=True)
            
            if st.button("Normalize", type="primary"):
                with st.spinner("Normalizing..."):
                    normalized, params = analyzer.normalize(
                        data, method=method, log_transform=log_transform
                    )
                
                st.session_state["normalized_data"] = normalized
                st.success("Normalization complete!")
                
                st.markdown("#### Normalized Data Preview")
                st.dataframe(normalized.head(10), use_container_width=True)
    
    with tab3:
        st.markdown("### PCA Analysis")
        
        if "normalized_data" not in st.session_state:
            st.warning("Please normalize data first.")
        else:
            data = st.session_state["normalized_data"]
            
            n_components = st.slider("Number of components", 2, 10, 2)
            
            if st.button("Run PCA", type="primary"):
                with st.spinner("Running PCA..."):
                    result, params = analyzer.pca(data, n_components=n_components)
                
                st.session_state["pca_result"] = result
                st.session_state["pca_params"] = params
                
                st.markdown("#### PCA Results")
                st.dataframe(result.head(), use_container_width=True)
                
                st.markdown("#### Explained Variance")
                for i, var in enumerate(params["explained_variance"]):
                    st.metric(f"PC{i+1}", f"{var*100:.1f}%")
                
                try:
                    fig_path = analyzer.plot_pca(
                        result,
                        params["explained_variance"][:2]
                    )
                    if fig_path:
                        st.image(str(fig_path))
                except Exception as e:
                    st.warning(f"Could not generate plot: {e}")
    
    with tab4:
        st.markdown("### Differential Expression")
        
        if "normalized_data" not in st.session_state:
            st.warning("Please normalize data first.")
        else:
            data = st.session_state["normalized_data"]
            
            st.markdown("#### Define Sample Groups")
            
            samples = list(data.columns)
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Group 1 (Treatment)**")
                group1 = st.multiselect("Select samples", samples, key="g1")
            
            with col2:
                st.markdown("**Group 2 (Control)**")
                group2 = st.multiselect("Select samples", samples, key="g2")
            
            fc_threshold = st.slider("Fold Change Threshold", 1.0, 5.0, 2.0)
            pval_threshold = st.slider("P-value Threshold", 0.001, 0.1, 0.05)
            
            if st.button("Run DE Analysis", type="primary"):
                if group1 and group2:
                    with st.spinner("Running DE analysis..."):
                        result, params = analyzer.differential_expression(
                            data,
                            {"treatment": group1, "control": group2},
                            fold_change_threshold=fc_threshold,
                            pvalue_threshold=pval_threshold,
                        )
                    
                    st.session_state["de_result"] = result
                    
                    st.markdown("#### Results Summary")
                    st.metric("Significant Genes", params["significant_genes"])
                    
                    st.markdown("#### Top Differentially Expressed Genes")
                    top_genes = result.nlargest(20, "log2_fold_change")
                    st.dataframe(top_genes, use_container_width=True)
                    
                    try:
                        fig_path = analyzer.plot_volcano(result)
                        if fig_path:
                            st.image(str(fig_path))
                    except Exception as e:
                        st.warning(f"Could not generate volcano plot: {e}")
                else:
                    st.error("Please select samples for both groups.")


def render_pipelines():
    st.markdown('<h1 class="subheader">📋 Pipeline Manager</h1>', unsafe_allow_html=True)
    
    st.markdown("### Available Pipelines")
    st.info("Pipeline management requires Snakemake or Nextflow to be installed.")
    
    st.markdown("#### Create New Pipeline")
    
    pipeline_name = st.text_input("Pipeline Name")
    pipeline_type = st.selectbox("Type", ["Snakemake", "Nextflow"])
    
    if pipeline_name:
        st.markdown("#### Configuration")
        
        config = {}
        config["data_dir"] = st.text_input("Data Directory", "data/raw")
        config["output_dir"] = st.text_input("Output Directory", "results")
        config["samples"] = st.text_area("Samples (comma-separated)", "sample1,sample2").split(",")
        
        if st.button("Create Pipeline", type="primary"):
            st.success(f"Pipeline '{pipeline_name}' configuration ready!")


def render_settings():
    st.markdown('<h1 class="subheader">⚙️ Settings</h1>', unsafe_allow_html=True)
    
    logger = get_logger()
    
    tab1, tab2, tab3 = st.tabs(["Environment", "Logs", "Export"])
    
    with tab1:
        st.markdown("### Environment Snapshot")
        
        if st.button("Capture Current Environment", type="primary"):
            snapshot = logger.capture_environment()
            st.success(f"Captured environment: {snapshot.hash}")
            st.json({
                "python_version": snapshot.python_version,
                "platform": snapshot.platform_info,
                "total_packages": len(snapshot.packages),
            })
        
        last_snapshot = logger.load_last_snapshot()
        if last_snapshot:
            st.markdown("#### Last Snapshot")
            st.json({
                "timestamp": last_snapshot.timestamp,
                "python_version": last_snapshot.python_version,
                "hash": last_snapshot.hash,
            })
    
    with tab2:
        st.markdown("### Command History")
        
        commands = logger.commands[-20:]
        if commands:
            for cmd in reversed(commands):
                with st.expander(f"{cmd.timestamp}: {cmd.command[:50]}..."):
                    st.json(asdict(cmd))
        else:
            st.info("No commands logged yet.")
    
    with tab3:
        st.markdown("### Export")
        
        export_format = st.selectbox("Format", ["json", "txt"])
        
        if st.button("Export Logs", type="primary"):
            output_path = logger.export_log("biopilot/logs/export.json", format=export_format)
            st.success(f"Exported to {output_path}")


def main():
    page = render_sidebar()
    
    if page == "🏠 Home":
        render_home()
    elif page == "🔍 Data Search":
        render_search()
    elif page == "📊 Sample Manager":
        render_samples()
    elif page == "📈 Analysis":
        render_analysis()
    elif page == "📋 Pipelines":
        render_pipelines()
    elif page == "⚙️ Settings":
        render_settings()


if __name__ == "__main__":
    main()
