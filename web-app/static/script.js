document.getElementById('analysisForm').addEventListener('submit', async function (e) {
    e.preventDefault();

    const submitBtn = document.getElementById('submitBtn');
    const originalBtnText = submitBtn.innerText;
    submitBtn.innerText = 'Analyzing...';
    submitBtn.disabled = true;

    // Gather data from form
    // Note: Inputs are strings by default, we need integers for the API
    const formData = {
        age: parseInt(document.getElementById('age').value),
        gender: parseInt(document.getElementById('gender').value),
        glucose: parseInt(document.getElementById('glucose').value),
        cholesterol: parseInt(document.getElementById('cholesterol').value),
        hdl: parseInt(document.getElementById('hdl').value),
        ldl: parseInt(document.getElementById('ldl').value),
        tch: parseInt(document.getElementById('tch').value),
        bmi: parseInt(document.getElementById('bmi').value),
        smoker: parseInt(document.getElementById('smoker').value),
        alcohol: parseInt(document.getElementById('alcohol').value)
    };

    try {
        const response = await fetch('/patient-data/predict', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(formData)
        });

        if (!response.ok) {
            throw new Error(`API Error: ${response.statusText}`);
        }

        const data = await response.json();
        displayResult(data.result);

    } catch (error) {
        console.error('Error:', error);
        alert('An error occurred while fetching the prediction. Please try again.');
    } finally {
        submitBtn.innerText = originalBtnText;
        submitBtn.disabled = false;
    }
});

function displayResult(score) {
    const container = document.getElementById('resultContainer');
    const scoreSpan = document.getElementById('riskScore');
    const bar = document.getElementById('riskBar');

    container.classList.remove('hidden');
    scoreSpan.innerText = score;

    // Animate the bar
    setTimeout(() => {
        bar.style.width = `${score}%`;

        // Dynamic color based on risk
        if (score < 30) {
            bar.style.backgroundColor = 'var(--success-color)';
            scoreSpan.style.color = 'var(--success-color)';
        } else if (score < 70) {
            bar.style.backgroundColor = 'var(--warning-color)';
            scoreSpan.style.color = 'var(--warning-color)';
        } else {
            bar.style.backgroundColor = 'var(--danger-color)';
            scoreSpan.style.color = 'var(--danger-color)';
        }
    }, 100);

    // Scroll view to result
    container.scrollIntoView({ behavior: 'smooth' });
}
