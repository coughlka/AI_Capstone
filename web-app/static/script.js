/**
 * CRC Biomarker Evidence Browser - Frontend Logic
 */

// State
let currentPage = 1;
let totalPages = 1;
let currentSort = 'final_score';
let currentSortOrder = 'desc';

// DOM Elements
const resultsBody = document.getElementById('resultsBody');
const pageInfo = document.getElementById('pageInfo');
const prevBtn = document.getElementById('prevPage');
const nextBtn = document.getElementById('nextPage');
const searchInput = document.getElementById('search');
const directionSelect = document.getElementById('direction');
const minScoreSlider = document.getElementById('minScore');
const minScoreValue = document.getElementById('minScoreValue');
const applyBtn = document.getElementById('applyFilters');
const resetBtn = document.getElementById('resetFilters');
const modal = document.getElementById('modal');
const modalClose = document.querySelector('.modal-close');
const modalTitle = document.getElementById('modalTitle');
const modalBody = document.getElementById('modalBody');

// Initialize
document.addEventListener('DOMContentLoaded', () => {
    loadStats();
    loadCandidates();
    setupEventListeners();
});

function setupEventListeners() {
    // Pagination
    prevBtn.addEventListener('click', () => {
        if (currentPage > 1) {
            currentPage--;
            loadCandidates();
        }
    });

    nextBtn.addEventListener('click', () => {
        if (currentPage < totalPages) {
            currentPage++;
            loadCandidates();
        }
    });

    // Filters
    applyBtn.addEventListener('click', () => {
        currentPage = 1;
        loadCandidates();
    });

    resetBtn.addEventListener('click', () => {
        searchInput.value = '';
        directionSelect.value = '';
        minScoreSlider.value = 0;
        minScoreValue.textContent = '0';
        currentPage = 1;
        currentSort = 'final_score';
        currentSortOrder = 'desc';
        loadCandidates();
    });

    // Min score slider
    minScoreSlider.addEventListener('input', () => {
        minScoreValue.textContent = minScoreSlider.value;
    });

    // Search on Enter
    searchInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
            currentPage = 1;
            loadCandidates();
        }
    });

    // Sortable headers
    document.querySelectorAll('th[data-sort]').forEach(th => {
        th.style.cursor = 'pointer';
        th.addEventListener('click', () => {
            const sortField = th.dataset.sort;
            if (currentSort === sortField) {
                currentSortOrder = currentSortOrder === 'asc' ? 'desc' : 'asc';
            } else {
                currentSort = sortField;
                currentSortOrder = 'desc';
            }
            currentPage = 1;
            loadCandidates();
        });
    });

    // Modal
    modalClose.addEventListener('click', closeModal);
    modal.addEventListener('click', (e) => {
        if (e.target === modal) closeModal();
    });
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape') closeModal();
    });
}

async function loadStats() {
    try {
        const response = await fetch('/api/stats');
        if (!response.ok) throw new Error('Failed to load stats');
        const stats = await response.json();

        document.getElementById('totalGenes').textContent =
            stats.total_genes?.toLocaleString() || '--';
        document.getElementById('significantGenes').textContent =
            stats.significant_genes?.toLocaleString() || '--';
        document.getElementById('scoreRange').textContent =
            stats.score_range ? `${stats.score_range.min.toFixed(1)} - ${stats.score_range.max.toFixed(1)}` : '--';
    } catch (error) {
        console.error('Error loading stats:', error);
    }
}

async function loadCandidates() {
    resultsBody.innerHTML = '<tr><td colspan="9" class="loading">Loading...</td></tr>';

    // Build query params
    const params = new URLSearchParams({
        page: currentPage,
        per_page: 50,
        sort_by: currentSort,
        sort_order: currentSortOrder
    });

    const search = searchInput.value.trim();
    if (search) params.append('search', search);

    const direction = directionSelect.value;
    if (direction) params.append('direction', direction);

    const minScore = parseInt(minScoreSlider.value);
    if (minScore > 0) params.append('min_score', minScore);

    try {
        const response = await fetch(`/api/candidates?${params}`);
        if (!response.ok) throw new Error('Failed to load candidates');
        const data = await response.json();

        renderTable(data.candidates);
        updatePagination(data.page, data.total_pages, data.total);
    } catch (error) {
        console.error('Error loading candidates:', error);
        resultsBody.innerHTML = `<tr><td colspan="9" class="error">Error loading data. Make sure the pipeline has been run.</td></tr>`;
    }
}

function renderTable(candidates) {
    if (!candidates || candidates.length === 0) {
        resultsBody.innerHTML = '<tr><td colspan="9" class="empty">No results found</td></tr>';
        return;
    }

    resultsBody.innerHTML = candidates.map(c => `
        <tr>
            <td class="gene-id" title="${c.gene}">${truncate(c.gene, 20)}</td>
            <td>${c.gene_symbol || '-'}</td>
            <td>
                <div class="score-cell">
                    <span class="score-value">${formatScore(c.final_score)}</span>
                    <div class="score-bar" style="width: ${c.final_score || 0}%"></div>
                </div>
            </td>
            <td>${formatScore(c.omics_score)}</td>
            <td>${formatScore(c.literature_score)}</td>
            <td>${formatScore(c.pathway_score)}</td>
            <td class="${getLog2fcClass(c.log2fc)}">${formatLog2fc(c.log2fc)}</td>
            <td>
                <span class="direction-badge ${c.direction || ''}">${c.direction || '-'}</span>
            </td>
            <td>
                <button class="btn-small" onclick="showGeneDetails('${c.gene}')">View</button>
            </td>
        </tr>
    `).join('');
}

function updatePagination(page, total, totalItems) {
    currentPage = page;
    totalPages = total || 1;
    pageInfo.textContent = `Page ${page} of ${totalPages} (${totalItems?.toLocaleString() || 0} genes)`;
    prevBtn.disabled = page <= 1;
    nextBtn.disabled = page >= totalPages;
}

async function showGeneDetails(geneId) {
    modalTitle.textContent = 'Loading...';
    modalBody.innerHTML = '<p>Loading gene details...</p>';
    modal.classList.remove('hidden');

    try {
        const response = await fetch(`/api/genes/${encodeURIComponent(geneId)}`);
        if (!response.ok) throw new Error('Failed to load gene details');
        const data = await response.json();

        modalTitle.textContent = data.gene_symbol ? `${data.gene_symbol} (${data.gene})` : data.gene;
        modalBody.innerHTML = renderGeneDetails(data);
    } catch (error) {
        console.error('Error loading gene details:', error);
        modalBody.innerHTML = '<p class="error">Error loading gene details</p>';
    }
}

function renderGeneDetails(data) {
    const scores = data.scores || {};
    const omics = data.omics_evidence || {};
    const pathway = data.pathway_evidence || {};

    return `
        <div class="detail-section">
            <h3>Evidence Scores</h3>
            <div class="scores-grid">
                <div class="score-item">
                    <span class="score-label">Final Score</span>
                    <span class="score-value-large">${formatScore(scores.final)}</span>
                </div>
                <div class="score-item">
                    <span class="score-label">Omics (45%)</span>
                    <span class="score-value-large">${formatScore(scores.omics)}</span>
                </div>
                <div class="score-item">
                    <span class="score-label">Literature (35%)</span>
                    <span class="score-value-large">${formatScore(scores.literature)}</span>
                </div>
                <div class="score-item">
                    <span class="score-label">Pathway (20%)</span>
                    <span class="score-value-large">${formatScore(scores.pathway)}</span>
                </div>
            </div>
        </div>

        <div class="detail-section">
            <h3>Omics Evidence</h3>
            ${omics.gene ? `
                <table class="detail-table">
                    <tr><td>Log2 Fold Change</td><td class="${getLog2fcClass(omics.log2fc)}">${formatLog2fc(omics.log2fc)}</td></tr>
                    <tr><td>Direction</td><td><span class="direction-badge ${omics.direction || ''}">${omics.direction || '-'}</span></td></tr>
                    <tr><td>FDR</td><td>${formatScientific(omics.fdr)}</td></tr>
                    <tr><td>P-value</td><td>${formatScientific(omics.p_value)}</td></tr>
                    <tr><td>Tumor Mean (log2)</td><td>${omics.tumor_mean?.toFixed(2) || '-'}</td></tr>
                    <tr><td>Normal Mean (log2)</td><td>${omics.normal_mean?.toFixed(2) || '-'}</td></tr>
                </table>
            ` : '<p>No omics evidence available</p>'}
        </div>

        <div class="detail-section">
            <h3>Pathway Evidence</h3>
            ${pathway.gene ? `
                <table class="detail-table">
                    <tr><td>Pathway Count</td><td>${pathway.pathway_count || 0}</td></tr>
                    <tr><td>Top Pathways</td><td>${pathway.top_pathways || '-'}</td></tr>
                </table>
            ` : '<p>No pathway evidence available</p>'}
        </div>
    `;
}

function closeModal() {
    modal.classList.add('hidden');
}

// Utility functions
function formatScore(value) {
    if (value === null || value === undefined) return '-';
    return value.toFixed(2);
}

function formatLog2fc(value) {
    if (value === null || value === undefined) return '-';
    const sign = value >= 0 ? '+' : '';
    return `${sign}${value.toFixed(2)}`;
}

function formatScientific(value) {
    if (value === null || value === undefined) return '-';
    if (value === 0) return '0';
    return value.toExponential(2);
}

function getLog2fcClass(value) {
    if (value === null || value === undefined) return '';
    return value >= 0 ? 'positive' : 'negative';
}

function truncate(str, maxLen) {
    if (!str) return '-';
    return str.length > maxLen ? str.substring(0, maxLen) + '...' : str;
}
