"""Data analyzer for normalization, PCA, and differential expression."""

import json
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from datetime import datetime
from dataclasses import dataclass, asdict
import warnings

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    from sklearn.manifold import TSNE
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    warnings.warn("scikit-learn not available. Some features will be disabled.")

try:
    from scipy import stats
    from scipy.cluster.hierarchy import linkage, dendrogram
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    warnings.warn("scipy not available. Some features will be disabled.")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    warnings.warn("matplotlib/seaborn not available. Plotting disabled.")


@dataclass
class AnalysisResult:
    analysis_id: str
    analysis_type: str
    parameters: Dict[str, Any]
    results: Dict[str, Any]
    created_at: str
    input_file: Optional[str] = None
    output_file: Optional[str] = None


class DataAnalyzer:
    """Perform normalization, PCA, and differential expression analysis."""
    
    def __init__(
        self,
        results_dir: str = "results",
        logs_dir: str = "logs",
    ):
        self.results_dir = Path(results_dir)
        self.figures_dir = self.results_dir / "figures"
        self.tables_dir = self.results_dir / "tables"
        self.logs_dir = Path(logs_dir)
        
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        
        self.analyses_file = self.logs_dir / "analyses.json"
        self.analyses: List[AnalysisResult] = self._load_analyses()
    
    def _load_analyses(self) -> List[AnalysisResult]:
        if not self.analyses_file.exists():
            return []
        with open(self.analyses_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        return [AnalysisResult(**a) for a in data]
    
    def _save_analyses(self):
        with open(self.analyses_file, "w", encoding="utf-8") as f:
            json.dump([asdict(a) for a in self.analyses], f, indent=2, default=str)
    
    def _generate_analysis_id(self) -> str:
        return f"analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    def load_expression_matrix(
        self,
        file_path: str,
        sep: str = "\t",
        gene_col: str = "gene",
        sample_cols: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        logger.info(f"Loading expression matrix from {file_path}")
        
        df = pd.read_csv(file_path, sep=sep)
        
        if gene_col in df.columns:
            df = df.set_index(gene_col)
        
        if sample_cols:
            df = df[sample_cols]
        
        df = df.select_dtypes(include=[np.number])
        
        logger.info(f"Loaded matrix: {df.shape[0]} genes x {df.shape[1]} samples")
        return df
    
    def normalize(
        self,
        data: pd.DataFrame,
        method: str = "cpm",
        log_transform: bool = True,
        **kwargs,
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Normalize expression data using various methods."""
        logger.info(f"Normalizing data using {method} method")
        
        params = {"method": method, "log_transform": log_transform, **kwargs}
        
        if method == "cpm":
            normalized = self._normalize_cpm(data, **kwargs)
        elif method == "rpkm":
            normalized = self._normalize_rpkm(data, **kwargs)
        elif method == "tpm":
            normalized = self._normalize_tpm(data, **kwargs)
        elif method == "zscore":
            normalized = self._normalize_zscore(data, **kwargs)
        elif method == "quantile":
            normalized = self._normalize_quantile(data, **kwargs)
        elif method == "minmax":
            normalized = self._normalize_minmax(data, **kwargs)
        elif method == "robust":
            normalized = self._normalize_robust(data, **kwargs)
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        if log_transform and method in ["cpm", "rpkm", "tpm"]:
            normalized = np.log2(normalized + 1)
            params["log2_transformed"] = True
        
        return normalized, params
    
    def _normalize_cpm(self, data: pd.DataFrame, **kwargs) -> pd.DataFrame:
        library_sizes = data.sum(axis=0)
        scaling_factor = kwargs.get("scaling_factor", 1e6)
        normalized = (data / library_sizes) * scaling_factor
        return normalized
    
    def _normalize_rpkm(
        self,
        data: pd.DataFrame,
        gene_lengths: Optional[Dict[str, float]] = None,
        **kwargs,
    ) -> pd.DataFrame:
        if gene_lengths is None:
            warnings.warn("No gene lengths provided. Using default normalization.")
            return self._normalize_cpm(data, **kwargs)
        
        lengths = pd.Series(gene_lengths)
        lengths = lengths.reindex(data.index)
        
        library_sizes = data.sum(axis=0)
        normalized = (data / library_sizes) * 1e6
        normalized = normalized.div(lengths, axis=0) * 1e3
        
        return normalized
    
    def _normalize_tpm(
        self,
        data: pd.DataFrame,
        gene_lengths: Optional[Dict[str, float]] = None,
        **kwargs,
    ) -> pd.DataFrame:
        if gene_lengths is None:
            warnings.warn("No gene lengths provided. Using CPM normalization.")
            return self._normalize_cpm(data, **kwargs)
        
        lengths = pd.Series(gene_lengths)
        lengths = lengths.reindex(data.index)
        
        rpk = data.div(lengths, axis=0) * 1e3
        rpk_sum = rpk.sum(axis=0)
        tpm = (rpk / rpk_sum) * 1e6
        
        return tpm
    
    def _normalize_zscore(self, data: pd.DataFrame, **kwargs) -> pd.DataFrame:
        if not SKLEARN_AVAILABLE:
            return (data - data.mean()) / data.std()
        
        scaler = StandardScaler(**kwargs)
        normalized = pd.DataFrame(
            scaler.fit_transform(data.T).T,
            index=data.index,
            columns=data.columns,
        )
        return normalized
    
    def _normalize_quantile(self, data: pd.DataFrame, **kwargs) -> pd.DataFrame:
        sorted_data = np.sort(data.values, axis=0)
        rank_mean = sorted_data.mean(axis=1)
        
        normalized = pd.DataFrame(index=data.index, columns=data.columns)
        for col in data.columns:
            ranks = data[col].rank(method="average").values.astype(int) - 1
            normalized[col] = rank_mean[ranks]
        
        return normalized
    
    def _normalize_minmax(self, data: pd.DataFrame, **kwargs) -> pd.DataFrame:
        if not SKLEARN_AVAILABLE:
            return (data - data.min()) / (data.max() - data.min())
        
        scaler = MinMaxScaler(**kwargs)
        normalized = pd.DataFrame(
            scaler.fit_transform(data.T).T,
            index=data.index,
            columns=data.columns,
        )
        return normalized
    
    def _normalize_robust(self, data: pd.DataFrame, **kwargs) -> pd.DataFrame:
        if not SKLEARN_AVAILABLE:
            median = data.median()
            q1 = data.quantile(0.25)
            q3 = data.quantile(0.75)
            iqr = q3 - q1
            return (data - median) / iqr
        
        scaler = RobustScaler(**kwargs)
        normalized = pd.DataFrame(
            scaler.fit_transform(data.T).T,
            index=data.index,
            columns=data.columns,
        )
        return normalized
    
    def pca(
        self,
        data: pd.DataFrame,
        n_components: int = 2,
        scale: bool = True,
        metadata: Optional[pd.DataFrame] = None,
        color_by: Optional[str] = None,
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Perform PCA analysis."""
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn required for PCA")
        
        logger.info(f"Performing PCA with {n_components} components")
        
        if scale:
            data_scaled = StandardScaler().fit_transform(data.T)
        else:
            data_scaled = data.T.values
        
        pca = PCA(n_components=n_components)
        components = pca.fit_transform(data_scaled)
        
        result_df = pd.DataFrame(
            components,
            index=data.columns,
            columns=[f"PC{i+1}" for i in range(n_components)],
        )
        
        if metadata is not None:
            result_df = result_df.join(metadata)
        
        params = {
            "n_components": n_components,
            "scale": scale,
            "explained_variance": pca.explained_variance_ratio_.tolist(),
            "total_variance_explained": sum(pca.explained_variance_ratio_),
        }
        
        return result_df, params
    
    def tsne(
        self,
        data: pd.DataFrame,
        n_components: int = 2,
        perplexity: float = 30.0,
        n_iter: int = 1000,
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Perform t-SNE analysis."""
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn required for t-SNE")
        
        logger.info("Performing t-SNE analysis")
        
        data_scaled = StandardScaler().fit_transform(data.T)
        
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            n_iter=n_iter,
            random_state=42,
        )
        embedding = tsne.fit_transform(data_scaled)
        
        result_df = pd.DataFrame(
            embedding,
            index=data.columns,
            columns=[f"tSNE{i+1}" for i in range(n_components)],
        )
        
        params = {
            "n_components": n_components,
            "perplexity": perplexity,
            "n_iter": n_iter,
        }
        
        return result_df, params
    
    def hierarchical_clustering(
        self,
        data: pd.DataFrame,
        method: str = "average",
        metric: str = "correlation",
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Perform hierarchical clustering."""
        if not SCIPY_AVAILABLE:
            raise ImportError("scipy required for hierarchical clustering")
        
        logger.info(f"Performing hierarchical clustering ({method}, {metric})")
        
        linkage_matrix = linkage(data.T, method=method, metric=metric)
        
        params = {"method": method, "metric": metric}
        
        return linkage_matrix, params
    
    def differential_expression(
        self,
        data: pd.DataFrame,
        sample_groups: Dict[str, List[str]],
        method: str = "ttest",
        fold_change_threshold: float = 2.0,
        pvalue_threshold: float = 0.05,
        multiple_testing: str = "fdr_bh",
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Perform differential expression analysis."""
        if not SCIPY_AVAILABLE:
            raise ImportError("scipy required for differential expression")
        
        logger.info(f"Performing differential expression: {method}")
        
        group_names = list(sample_groups.keys())
        if len(group_names) != 2:
            raise ValueError("Exactly two groups required for differential expression")
        
        group1 = sample_groups[group_names[0]]
        group2 = sample_groups[group_names[1]]
        
        group1_data = data[group1]
        group2_data = data[group2]
        
        results = []
        
        for gene in data.index:
            g1_vals = group1_data.loc[gene].values
            g2_vals = group2_data.loc[gene].values
            
            if method == "ttest":
                stat, pvalue = stats.ttest_ind(g1_vals, g2_vals)
            elif method == "ranksum":
                stat, pvalue = stats.ranksums(g1_vals, g2_vals)
            elif method == "mannwhitney":
                stat, pvalue = stats.mannwhitneyu(g1_vals, g2_vals, alternative="two-sided")
            else:
                raise ValueError(f"Unknown method: {method}")
            
            mean_g1 = np.mean(g1_vals)
            mean_g2 = np.mean(g2_vals)
            
            log2_fc = np.log2((mean_g1 + 1e-10) / (mean_g2 + 1e-10))
            
            results.append({
                "gene": gene,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
                "log2_fold_change": log2_fc,
                "statistic": stat,
                "pvalue": pvalue,
            })
        
        result_df = pd.DataFrame(results)
        result_df = result_df.set_index("gene")
        
        if multiple_testing and SCIPY_AVAILABLE:
            from statsmodels.stats.multitest import multipletests
            rejected, pvals_corrected, _, _ = multipletests(
                result_df["pvalue"].values,
                method=multiple_testing,
            )
            result_df["adj_pvalue"] = pvals_corrected
            result_df["significant"] = (
                (result_df["adj_pvalue"] < pvalue_threshold) &
                (np.abs(result_df["log2_fold_change"]) >= np.log2(fold_change_threshold))
            )
        else:
            result_df["significant"] = (
                (result_df["pvalue"] < pvalue_threshold) &
                (np.abs(result_df["log2_fold_change"]) >= np.log2(fold_change_threshold))
            )
        
        params = {
            "method": method,
            "groups": sample_groups,
            "fold_change_threshold": fold_change_threshold,
            "pvalue_threshold": pvalue_threshold,
            "multiple_testing": multiple_testing,
            "total_genes": len(result_df),
            "significant_genes": result_df["significant"].sum(),
        }
        
        return result_df, params
    
    def plot_pca(
        self,
        pca_result: pd.DataFrame,
        explained_variance: List[float],
        color_by: Optional[str] = None,
        output_file: Optional[str] = None,
    ) -> Optional[Path]:
        """Create PCA scatter plot."""
        if not PLOTTING_AVAILABLE:
            warnings.warn("Plotting not available")
            return None
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        if color_by and color_by in pca_result.columns:
            groups = pca_result[color_by].unique()
            for group in groups:
                mask = pca_result[color_by] == group
                ax.scatter(
                    pca_result.loc[mask, "PC1"],
                    pca_result.loc[mask, "PC2"],
                    label=group,
                    alpha=0.7,
                )
            ax.legend()
        else:
            ax.scatter(pca_result["PC1"], pca_result["PC2"], alpha=0.7)
        
        ax.set_xlabel(f"PC1 ({explained_variance[0]*100:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({explained_variance[1]*100:.1f}% variance)")
        ax.set_title("PCA Analysis")
        
        output_path = Path(output_file) if output_file else self.figures_dir / "pca_plot.png"
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        logger.info(f"Saved PCA plot to {output_path}")
        return output_path
    
    def plot_volcano(
        self,
        de_result: pd.DataFrame,
        log2fc_col: str = "log2_fold_change",
        pval_col: str = "adj_pvalue",
        output_file: Optional[str] = None,
    ) -> Optional[Path]:
        """Create volcano plot."""
        if not PLOTTING_AVAILABLE:
            warnings.warn("Plotting not available")
            return None
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        log2fc = de_result[log2fc_col]
        neg_log_pval = -np.log10(de_result[pval_col])
        
        colors = np.where(
            de_result["significant"],
            np.where(log2fc > 0, "red", "blue"),
            "gray",
        )
        
        ax.scatter(log2fc, neg_log_pval, c=colors, alpha=0.6, s=20)
        
        ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.5)
        ax.axvline(-1, color="black", linestyle="--", linewidth=0.5)
        ax.axvline(1, color="black", linestyle="--", linewidth=0.5)
        
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10 Adjusted P-value")
        ax.set_title("Volcano Plot")
        
        output_path = Path(output_file) if output_file else self.figures_dir / "volcano_plot.png"
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        logger.info(f"Saved volcano plot to {output_path}")
        return output_path
    
    def plot_heatmap(
        self,
        data: pd.DataFrame,
        metadata: Optional[pd.DataFrame] = None,
        cluster_genes: bool = True,
        cluster_samples: bool = True,
        output_file: Optional[str] = None,
        top_genes: int = 100,
    ) -> Optional[Path]:
        """Create expression heatmap."""
        if not PLOTTING_AVAILABLE:
            warnings.warn("Plotting not available")
            return None
        
        if data.shape[0] > top_genes:
            variances = data.var(axis=1)
            top_idx = variances.nlargest(top_genes).index
            data = data.loc[top_idx]
        
        fig_size = max(8, data.shape[1] * 0.5)
        fig, ax = plt.subplots(figsize=(fig_size, 12))
        
        cluster_kws = {}
        if cluster_genes:
            cluster_kws["row_cluster"] = True
        if cluster_samples:
            cluster_kws["col_cluster"] = True
        
        if metadata is not None:
            sample_colors = metadata.apply(lambda x: pd.Categorical(x).codes)
            cluster_kws["col_colors"] = sample_colors.values
        
        sns.clustermap(
            data,
            cmap="RdBu_r",
            center=0,
            figsize=(fig_size, 12),
            **cluster_kws,
        )
        
        output_path = Path(output_file) if output_file else self.figures_dir / "heatmap.png"
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        
        logger.info(f"Saved heatmap to {output_path}")
        return output_path
    
    def save_results(
        self,
        result: pd.DataFrame,
        analysis_type: str,
        params: Dict[str, Any],
        input_file: Optional[str] = None,
    ) -> Tuple[Path, AnalysisResult]:
        analysis_id = self._generate_analysis_id()
        output_file = self.tables_dir / f"{analysis_id}.csv"
        
        result.to_csv(output_file)
        
        analysis = AnalysisResult(
            analysis_id=analysis_id,
            analysis_type=analysis_type,
            parameters=params,
            results={"output_file": str(output_file)},
            created_at=datetime.now().isoformat(),
            input_file=input_file,
            output_file=str(output_file),
        )
        
        self.analyses.append(analysis)
        self._save_analyses()
        
        logger.info(f"Saved results to {output_file}")
        return output_file, analysis
    
    def get_analysis(self, analysis_id: str) -> Optional[AnalysisResult]:
        for analysis in self.analyses:
            if analysis.analysis_id == analysis_id:
                return analysis
        return None
    
    def list_analyses(self, analysis_type: Optional[str] = None) -> List[AnalysisResult]:
        if analysis_type:
            return [a for a in self.analyses if a.analysis_type == analysis_type]
        return self.analyses
